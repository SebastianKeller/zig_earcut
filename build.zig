const Builder = @import("std").build.Builder;

pub fn build(b: *Builder) void {
    const mode = b.standardReleaseOptions();
    const lib = b.addStaticLibrary("zig_earcut", "src/main.zig");
    lib.setBuildMode(mode);
    lib.install();

    var tests = b.addTest("test/test.zig");
    tests.addPackagePath("earcut", "src/main.zig");
    tests.setBuildMode(mode);
    tests.linkLibC();

    const test_step = b.step("test", "Run library tests");
    test_step.dependOn(&tests.step);
}
