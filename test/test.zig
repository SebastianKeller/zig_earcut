const std = @import("std");
const Earcut = @import("earcut").Earcut;

test "flattening" {
    var rings = try std.testing.allocator.alloc([][2]f64, 1);
    var points = try std.testing.allocator.alloc([2]f64, 4);
    points[0] = .{ 10, 0 };
    points[1] = .{ 0, 50 };
    points[2] = .{ 60, 60 };
    points[3] = .{ 70, 10 };
    rings[0] = points;

    defer std.testing.allocator.free(points);
    defer std.testing.allocator.free(rings);

    var expected = [_]f64{ 10, 0, 0, 50, 60, 60, 70, 10 };

    var flattened = try Earcut.flatten(2, rings, std.testing.allocator);
    std.testing.expectEqualSlices(f64, &expected, flattened.vertices);
    flattened.deinit();
}

test "indices-2d" {
    var earcut = Earcut.init(std.testing.allocator);
    defer earcut.deinit();

    var points = [_]f64{
        10, 0,
        0,  50,
        60, 60,
        70, 10,
    };
    var result = try earcut.earcut(&points, null, 2);
    std.testing.expectEqualSlices(usize, result, &[_]usize{ 1, 0, 3, 3, 2, 1 });
}

test "indices-3d" {
    var earcut = Earcut.init(std.testing.allocator);
    defer earcut.deinit();

    var points = [_]f64{
        10, 0,  0,
        0,  50, 0,
        60, 60, 0,
        70, 10, 0,
    };
    var result = try earcut.earcut(&points, null, 3);
    std.testing.expectEqualSlices(usize, result, &[_]usize{ 1, 0, 3, 3, 2, 1 });
}
