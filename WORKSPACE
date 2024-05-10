load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

http_archive(
  name = "com_google_googletest",
  urls = ["https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip"],
  #urls = ["https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip"],
  strip_prefix = "googletest-1.14.0",
)

#http_archive(
#    name = "random_events",
#    urls = ["https://github.com/tomsch420/random-events-lib/archive/refs/tags/0.0.1.zip"],
#    build_file = "//third_party:random_events.BUILD.bazel",
#)

new_local_repository(
    name = "random_events",
    path = "/home/tom_sch/random-events-lib",
    build_file = "//third_party:random_events.BUILD.bazel",
)