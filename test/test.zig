const std = @import("std");
const Earcut = @import("earcut").Earcut;

test "flattening" {
    inline for ([_]type{ f32, f64, f128 }) |Scalar| {
        var rings = try std.testing.allocator.alloc([][2]Scalar, 1);
        var points = try std.testing.allocator.alloc([2]Scalar, 4);
        points[0] = .{ 10, 0 };
        points[1] = .{ 0, 50 };
        points[2] = .{ 60, 60 };
        points[3] = .{ 70, 10 };
        rings[0] = points;

        defer std.testing.allocator.free(points);
        defer std.testing.allocator.free(rings);

        var expected = [_]Scalar{ 10, 0, 0, 50, 60, 60, 70, 10 };

        var flattened = try Earcut(Scalar).flatten(2, rings, std.testing.allocator);
        std.testing.expectEqualSlices(Scalar, &expected, flattened.vertices);
        flattened.deinit();
    }
}

test "indices-2d" {
    inline for ([_]type{ f32, f64, f128 }) |Scalar| {
        var earcut = Earcut(Scalar).init(std.testing.allocator);
        defer earcut.deinit();

        var points = [_]Scalar{
            10, 0,
            0,  50,
            60, 60,
            70, 10,
        };
        var result = try earcut.earcut(&points, null, 2);
        std.testing.expectEqualSlices(usize, result, &[_]usize{ 1, 0, 3, 3, 2, 1 });
    }
}

test "indices-3d" {
    inline for ([_]type{ f32, f64, f128 }) |Scalar| {
        var earcut = Earcut(Scalar).init(std.testing.allocator);
        defer earcut.deinit();

        var points = [_]Scalar{
            10, 0,  0,
            0,  50, 0,
            60, 60, 0,
            70, 10, 0,
        };
        var result = try earcut.earcut(&points, null, 3);
        std.testing.expectEqualSlices(usize, result, &[_]usize{ 1, 0, 3, 3, 2, 1 });
    }
}
