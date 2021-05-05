const std = @import("std");

pub fn Earcut(comptime Scalar: type) type {
    const Scalar_min = switch (Scalar) {
        f16 => std.math.f16_min,
        f32 => std.math.f32_min,
        f64 => std.math.f64_min,
        f128 => std.math.f128_min,
        else => unreachable,
    };

    const Scalar_max = switch (Scalar) {
        f16 => std.math.f16_max,
        f32 => std.math.f32_max,
        f64 => std.math.f64_max,
        f128 => std.math.f128_max,
        else => unreachable,
    };

    const Node = struct {
        i: usize,
        x: Scalar,
        y: Scalar,
        z: Scalar = Scalar_min,
        steiner: bool = false,

        prev: *@This(),
        next: *@This(),
        prevZ: ?*@This() = null,
        nextZ: ?*@This() = null,
    };

    return struct {
        arena: std.heap.ArenaAllocator,

        const Self = @This();

        pub fn init(allocator: *std.mem.Allocator) Self {
            return .{
                .arena = std.heap.ArenaAllocator.init(allocator),
            };
        }

        pub fn deinit(self: Self) void {
            self.arena.deinit();
        }

        pub fn earcut(self: *Self, data: []Scalar, hole_indicies: ?[]usize, dim: usize) ![]usize {
            const hasHoles = hole_indicies != null and hole_indicies.?.len > 0;
            const outerLen = if (hasHoles) hole_indicies.?[0] * dim else data.len;

            const outerNode = try self.linkedList(data, 0, outerLen, dim, true);

            var triangles = std.ArrayList(usize).init(&self.arena.allocator);
            if (outerNode == null)
                return triangles.items;

            var minX: Scalar = 0.0;
            var minY: Scalar = 0.0;
            var maxX: Scalar = 0.0;
            var maxY: Scalar = 0.0;
            var size: Scalar = Scalar_min;

            var node = outerNode.?;
            if (hasHoles)
                node = try self.eliminateHoles(data, hole_indicies.?, node, dim);

            // if the shape is not too simple, we'll use z-order curve hash later;
            // calculate polygon bbox
            if (data.len > 80 * dim) {
                minX = data[0];
                maxX = data[0];

                minY = data[1];
                maxY = data[1];

                var i = dim;
                while (i < outerLen) : (i += dim) {
                    const x = data[i];
                    const y = data[i + 1];

                    if (x < minX) minX = x;
                    if (y < minY) minY = y;

                    if (x > maxX) maxX = x;
                    if (y > maxY) maxY = y;
                }

                // minX, minY and size are later used to transform coords into integers for z-order calculation
                size = std.math.max(maxX - minX, maxY - minY);
            }

            try self.earcutLinked(node, &triangles, dim, minX, minY, size, -1);
            return triangles.items;
        }

        pub const Flattened = struct {
            allocator: *std.mem.Allocator,
            vertices: []Scalar,
            holes: []usize,
            dimension: usize,

            pub fn deinit(self: Flattened) void {
                self.allocator.free(self.vertices);
                self.allocator.free(self.holes);
            }
        };

        pub fn flatten(comptime dim: usize, data: [][][dim]Scalar, allocator: *std.mem.Allocator) !Flattened {
            var vertices_count: usize = 0;
            var hole_count: usize = data.len - 1;
            for (data) |ring| {
                vertices_count += ring.len * dim;
            }

            var vertices = try allocator.alloc(Scalar, vertices_count);
            var hole_indexes = try allocator.alloc(usize, hole_count);

            var vert_idx: usize = 0;
            var hole_idx: usize = 0;
            var i: usize = 0;
            while (i < data.len) : (i += 1) {
                if (i > 0) {
                    hole_indexes[hole_idx] = vert_idx / dim;
                    hole_idx += 1;
                }

                var j: usize = 0;
                while (j < data[i].len) : (j += 1) {
                    var d: usize = 0;
                    while (d < dim) : (d += 1) {
                        vertices[vert_idx] = data[i][j][d];
                        vert_idx += 1;
                    }
                }
            }

            return Flattened{ .allocator = allocator, .vertices = vertices, .holes = hole_indexes, .dimension = dim };
        }

        fn earcutLinked(
            self: *Self,
            ear: ?*Node,
            triangles: *std.ArrayList(usize),
            dim: usize,
            minX: Scalar,
            minY: Scalar,
            size: Scalar,
            pass: isize,
        ) std.mem.Allocator.Error!void {
            if (ear == null)
                return;

            if (pass == -1 and size != Scalar_min)
                _ = indexCurve(ear.?, minX, minY, size);

            var e = ear.?;
            var stop = e;

            while (e.prev != e.next) {
                const prev = e.prev;
                const next = e.next;

                if (if (size != Scalar_min) isEarHashed(e, minX, minY, size) else isEar(e)) {
                    try triangles.append(prev.i / dim);
                    try triangles.append(e.i / dim);
                    try triangles.append(next.i / dim);

                    removeNode(e);

                    e = next.next;
                    stop = next.next;
                    continue;
                }

                e = next;

                // if we looped through the while remining polygon and can't find any more ears
                if (e == stop) {
                    if (pass == -1) {
                        try self.earcutLinked(filterPoints(e, null), triangles, dim, minX, minY, size, 1);
                    }
                    // if this didn't work, try curing all small self-intersections locally
                    else if (pass == 1) {
                        e = try self.cureLocalIntersections(e, triangles, dim);
                        try self.earcutLinked(e, triangles, dim, minX, minY, size, 2);
                    }
                    // as a last resort, try splitting the remainig polygon into two
                    else if (pass == 2) {
                        try self.splitEarcut(e, triangles, dim, minX, minY, size);
                    }
                    break;
                }
            }
        }

        fn splitEarcut(
            self: *Self,
            start: *Node,
            triangles: *std.ArrayList(usize),
            dim: usize,
            minX: Scalar,
            minY: Scalar,
            size: Scalar,
        ) !void {
            var a = start;
            while (true) {
                var b = a.next.next;
                while (b != a.prev) {
                    if (a.i != b.i and isValidDiagonal(a, b)) {
                        var c = try self.splitPolygon(a, b);

                        //filter colinear points around the cuts
                        var wrapper: ?*Node = a.next;
                        a = filterPoints(a, wrapper).?;

                        wrapper = c.next;
                        c = filterPoints(c, wrapper).?;

                        // run earcut on each half
                        try self.earcutLinked(a, triangles, dim, minX, minY, size, -1);
                        try self.earcutLinked(c, triangles, dim, minX, minY, size, -1);
                        return;
                    }
                    b = b.next;
                }

                a = a.next;
                if (a == start)
                    return;
            }
        }

        fn cureLocalIntersections(self: *Self, start: *Node, triangles: *std.ArrayList(usize), dim: usize) !*Node {
            var s = start;
            var p = s;
            while (true) {
                var a = p.prev;
                var b = p.next.next;

                if (!equals(a, b) and intersects(a, p, p.next, b) and locallyInside(a, b) and locallyInside(b, a)) {
                    try triangles.append(a.i / dim);
                    try triangles.append(p.i / dim);
                    try triangles.append(b.i / dim);

                    // remove two nodes involved
                    removeNode(p);
                    removeNode(p.next);

                    p = b;
                    s = b;
                }
                p = p.next;

                if (p == s)
                    return p;
            }
        }

        fn isEar(ear: *const Node) bool {
            var a = ear.prev;
            var b = ear;
            var c = ear.next;

            if (area(a, b, c) >= 0)
                return false; //reflex, can't be an ear

            // now make sure we don't have other points inside the potential ear
            //
            var p = ear.next.next;

            while (p != ear.prev) {
                if (pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) and area(p.prev, p, p.next) >= 0)
                    return false;
                p = p.next;
            }
            return true;
        }

        fn isEarHashed(ear: *const Node, minX: Scalar, minY: Scalar, size: Scalar) bool {
            var a = ear.prev;
            var b = ear;
            var c = ear.next;

            if (area(a, b, c) >= 0)
                return false; // reflex, can't be an ear

            // triangle bbox; min & max are calculated like this for speed
            var minTX: Scalar = 0;
            var minTY: Scalar = 0;
            var maxTX: Scalar = 0;
            var maxTY: Scalar = 0;

            if (a.x < b.x) {
                minTX = if (a.x < c.x) a.x else c.x;
            } else {
                minTX = if (b.x < c.x) b.x else c.x;
            }

            if (a.y < b.y) {
                minTY = if (a.y < c.y) a.y else c.y;
            } else {
                minTY = if (b.y < c.y) b.y else c.y;
            }

            if (a.x > b.x) {
                maxTX = if (a.x > c.x) a.x else c.x;
            } else {
                maxTX = if (b.x > c.x) b.x else c.x;
            }

            if (a.y > b.y) {
                maxTY = if (a.y > c.y) a.y else c.y;
            } else {
                maxTY = if (b.y > c.y) b.y else c.y;
            }

            // z-order range for the current triangle bbox;
            const minZ = zOrder(minTX, minTY, minX, minY, size);
            const maxZ = zOrder(maxTX, maxTY, minX, minY, size);

            // first look for points inside the triangle in increasing z-order
            var p = ear.nextZ;

            while (p != null and p.?.z <= maxZ) {
                if (p != ear.prev and p != ear.next and pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.?.x, p.?.y) and area(p.?.prev, p.?, p.?.next) >= 0)
                    return false;
                p = p.?.nextZ;
            }

            // then look for points in decreasing z-order
            p = ear.prevZ;

            while (p != null and p.?.z >= minZ) {
                if (p != ear.prev and p != ear.next and pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.?.x, p.?.y) and area(p.?.prev, p.?, p.?.next) >= 0)
                    return false;
                p = p.?.prevZ;
            }

            return true;
        }

        fn zOrder(x: Scalar, y: Scalar, minX: Scalar, minY: Scalar, size: Scalar) Scalar {
            // coords are transformed into non-negative 15-bit integer range
            var lx = @floatToInt(i32, 32767 * (x - minX) / size);
            var ly = @floatToInt(i32, 32767 * (y - minY) / size);

            lx = (lx | (lx << 8)) & 0x00FF00FF;
            lx = (lx | (lx << 4)) & 0x0F0F0F0F;
            lx = (lx | (lx << 2)) & 0x33333333;
            lx = (lx | (lx << 1)) & 0x55555555;

            ly = (ly | (ly << 8)) & 0x00FF00FF;
            ly = (ly | (ly << 4)) & 0x0F0F0F0F;
            ly = (ly | (ly << 2)) & 0x33333333;
            ly = (ly | (ly << 1)) & 0x55555555;

            return @intToFloat(Scalar, lx | (ly << 1));
        }

        fn indexCurve(start: *Node, minX: Scalar, minY: Scalar, size: Scalar) void {
            var p = start;
            while (true) {
                if (p.z == Scalar_min)
                    p.z = zOrder(p.x, p.y, minX, minY, size);
                p.prevZ = p.prev;
                p.nextZ = p.next;
                p = p.next;
                if (p == start)
                    break;
            }

            p.prevZ.?.nextZ = null;
            p.prevZ = null;

            _ = sortLinked(p);
        }

        fn sortLinked(list: *Node) *Node {
            var inSize: i32 = 1;

            var l: ?*Node = list;
            while (true) {
                var p = l;
                l = null;
                var tail: ?*Node = null;
                var numMerges: usize = 0;

                while (p != null) {
                    numMerges += 1;
                    var q = p;
                    var pSize: usize = 0;

                    var i: usize = 0;
                    while (i < inSize) : (i += 1) {
                        pSize += 1;
                        q = q.?.nextZ;
                        if (q == null)
                            break;
                    }

                    var qSize = inSize;

                    while (pSize > 0 or (qSize > 0 and q != null)) {
                        var e: ?*Node = p;
                        if (pSize == 0) {
                            e = q;
                            q = q.?.nextZ;
                            qSize -= 1;
                        } else if (qSize == 0 or q == null) {
                            e = p;
                            p = p.?.nextZ;
                            pSize -= 1;
                        } else if (p.?.z <= q.?.z) {
                            e = p;
                            p = p.?.nextZ;
                            pSize -= 1;
                        } else {
                            e = q;
                            q = q.?.nextZ;
                            qSize -= 1;
                        }

                        if (tail != null) {
                            tail.?.nextZ = e;
                        } else {
                            l = e;
                        }

                        e.?.prevZ = tail;
                        tail = e;
                    }

                    p = q;
                }

                tail.?.nextZ = null;
                inSize *= 2;

                if (numMerges <= 1)
                    break;
            }

            return l.?;
        }

        fn eliminateHoles(self: *Self, data: []Scalar, holeIndices: []usize, outerNode: *Node, dim: usize) !*Node {
            var queue = std.ArrayList(*Node).init(&self.arena.allocator);

            const len = holeIndices.len;
            var i: usize = 0;
            while (i < len) : (i += 1) {
                var start = holeIndices[i] * dim;
                var end = if (i < len - 1) holeIndices[i + 1] * dim else data.len;
                var list = try self.linkedList(data, start, end, dim, false);
                if (list) |v| {
                    if (v == v.next)
                        v.steiner = true;

                    try queue.append(getLeftmost(v));
                }
            }

            std.sort.sort(*Node, queue.items, {}, nodeCompare);

            var n = outerNode;
            for (queue.items) |node| {
                try self.eliminateHole(node, n);
                n = filterPoints(n, n.next).?;
            }
            return n;
        }

        fn nodeCompare(context: void, left: *Node, right: *Node) bool {
            if (left.x > right.x)
                return false;
            return true;
        }

        fn filterPoints(start: *Node, end: ?*Node) ?*Node {
            var e = end;
            if (end == null)
                e = start;

            var p = start;
            while (true) {
                var again = false;
                if (!p.steiner and equals(p, p.next) or area(p.prev, p, p.next) == 0) {
                    removeNode(p);
                    p = p.prev;
                    e = p.prev;
                    if (p == p.next)
                        return null;
                    again = true;
                } else {
                    p = p.next;
                }

                if (!(again or p != e))
                    break;
            }

            return e;
        }

        fn eliminateHole(self: *Self, hole: *Node, outerNode: *Node) !void {
            var node = findHoleBridge(hole, outerNode);
            if (node) |n| {
                var b = try self.splitPolygon(n, hole);
                _ = filterPoints(b, b.next);
            }
        }

        fn splitPolygon(self: *Self, a: *Node, b: *Node) !*Node {
            var a2 = try self.arena.allocator.create(Node);
            a2.* = .{ .i = a.i, .x = a.x, .y = a.y, .next = a2, .prev = a2 };

            var b2 = try self.arena.allocator.create(Node);
            b2.* = .{ .i = b.i, .x = b.x, .y = b.y, .next = b2, .prev = b2 };

            var an = a.next;
            var bp = b.prev;

            a.next = b;
            b.prev = a;

            a2.next = an;
            an.prev = a2;

            b2.next = a2;
            a2.prev = b2;

            bp.next = b2;
            b2.prev = bp;

            return b2;
        }

        fn findHoleBridge(hole: *Node, outerNode: *Node) ?*Node {
            var p = outerNode;
            var hx = hole.x;
            var hy = hole.y;
            var qx: Scalar = -Scalar_max;
            var m: ?*Node = null;

            // find a segment intersected by a ray from the hole's leftmost point to
            // the left;
            // segment's endpoint with lesser x will be potential connection point
            while (true) {
                if (hy <= p.y and hy >= p.next.y) {
                    var x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y);
                    if (x <= hx and x > qx) {
                        qx = x;
                        if (x == hx) {
                            if (hy == p.y)
                                return p;
                            if (hy == p.next.y)
                                return p.next;
                        }
                        m = if (p.x < p.next.x) p else p.next;
                    }
                }
                p = p.next;
                if (p == outerNode)
                    break;
            }

            if (m == null)
                return null;

            if (hx == qx)
                return m.?.prev; // hole touches outer segment; pick lower endpoint

            // look for points inside the triangle of hole point, segment
            // intersection and endpoint;
            // if there are no points found, we have a valid connection;
            // otherwise choose the point of the minimum angle with the ray as
            // connection point

            var stop = m;
            var mx = m.?.x;
            var my = m.?.y;
            var tanMin: Scalar = Scalar_max;

            p = m.?.next;

            while (p != stop) {
                if (hx >= p.x and p.x >= mx and pointInTriangle(if (hy < my) hx else qx, hy, mx, my, if (hy < my) qx else hx, hy, p.x, p.y)) {
                    const tan = std.math.fabs(hy - p.y) / (hx - p.x); // tangential

                    if ((tan < tanMin or (tan == tanMin and p.x > m.?.x)) and locallyInside(p, hole)) {
                        m = p;
                        tanMin = tan;
                    }
                }

                p = p.next;
            }

            return m;
        }

        fn linkedList(self: *Self, data: []Scalar, start: usize, end: usize, dim: usize, clockwise: bool) !?*Node {
            var last: ?*Node = null;
            if (clockwise == (signedArea(data, start, end, dim) > 0)) {
                var i = start;
                while (i < end) : (i += dim) {
                    last = try self.insertNode(i, data[i], data[i + 1], last);
                }
            } else {
                var i = @intCast(i64, (end - dim));
                while (i >= start) : (i -= @intCast(i64, dim)) {
                    const idx = @intCast(usize, i);
                    last = try self.insertNode(idx, data[idx], data[idx + 1], last);
                }
            }

            if (last) |v| {
                if (equals(v, v.next)) {
                    removeNode(v);
                    last = v.next;
                }
            }
            return last;
        }

        fn removeNode(p: *Node) void {
            p.next.prev = p.prev;
            p.prev.next = p.next;

            if (p.prevZ) |v| {
                v.nextZ = p.nextZ;
            }

            if (p.nextZ) |v| {
                v.prevZ = p.prevZ;
            }
        }

        fn insertNode(self: *Self, i: usize, x: Scalar, y: Scalar, last: ?*Node) !*Node {
            var p = try self.arena.allocator.create(Node);

            if (last) |v| {
                p.* = .{ .i = i, .x = x, .y = y, .next = v.next, .prev = v };
                v.next.prev = p;
                v.next = p;
            } else {
                p.* = .{ .i = i, .x = x, .y = y, .next = p, .prev = p };
            }
            return p;
        }

        fn signedArea(data: []Scalar, start: usize, end: usize, dim: usize) Scalar {
            var sum: Scalar = 0.0;
            var j = end - dim;
            var i = start;
            while (i < end) : (i += dim) {
                sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
                j = i;
            }
            return sum;
        }

        fn getLeftmost(start: *Node) *Node {
            var p = start;
            var leftmost = start;
            while (true) {
                if (p.x < leftmost.x)
                    leftmost = p;
                p = p.next;
                if (p == start)
                    break;
            }
            return leftmost;
        }

        fn equals(p1: *const Node, p2: *const Node) bool {
            return p1.x == p2.x and p1.y == p2.y;
        }

        fn area(p: *const Node, q: *const Node, r: *const Node) Scalar {
            return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
        }

        fn isValidDiagonal(a: *const Node, b: *const Node) bool {
            return a.next.i != b.i and a.prev.i != b.i and !intersectsPolygon(a, b) and locallyInside(a, b) and locallyInside(b, a) and middleInside(a, b);
        }

        fn middleInside(a: *const Node, b: *const Node) bool {
            var p = a;
            var inside = false;
            var px = (a.x + b.x) / 2;
            var py = (a.y + b.y) / 2;
            while (true) {
                if (((p.y > py) != (p.next.y > py)) and (px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x))
                    inside = !inside;
                p = p.next;
                if (p == a)
                    break;
            }

            return inside;
        }

        fn intersectsPolygon(a: *const Node, b: *const Node) bool {
            var p = a;
            while (true) {
                if (p.i != a.i and p.next.i != a.i and p.i != b.i and p.next.i != b.i and intersects(p, p.next, a, b))
                    return true;
                p = p.next;
                if (p == a)
                    break;
            }
            return false;
        }

        fn intersects(p1: *const Node, q1: *const Node, p2: *const Node, q2: *const Node) bool {
            if ((equals(p1, q1) and equals(p2, q2)) or (equals(p1, q2) and equals(p2, q1)))
                return true;
            return (area(p1, q1, p2) > 0) != (area(p1, q1, q2) > 0) and (area(p2, q2, p1) > 0) != (area(p2, q2, q1) > 0);
        }

        fn locallyInside(a: *const Node, b: *const Node) bool {
            return if (area(a.prev, a, a.next) < 0)
                area(a, b, a.next) >= 0 and area(a, a.prev, b) >= 0
            else
                area(a, b, a.prev) < 0 or area(a, a.next, b) < 0;
        }

        fn pointInTriangle(ax: Scalar, ay: Scalar, bx: Scalar, by: Scalar, cx: Scalar, cy: Scalar, px: Scalar, py: Scalar) bool {
            return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 and (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 and (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0;
        }
    };
}
