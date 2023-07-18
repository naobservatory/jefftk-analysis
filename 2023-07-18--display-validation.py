#!/usr/bin/env python3
import sys
import fcntl
import termios
import struct

COLS = struct.unpack(
  "hh", fcntl.ioctl(0, termios.TIOCGWINSZ, "1234")
)[1]

COLOR_END = '\x1b[0m'
COLOR_RED = '\x1b[1;31m'

fname, = sys.argv[1:]

with open(fname) as inf:
  lines = [line[:-1] for line in inf]

longest_line_length = max(len(line) for line in lines)

skip_start = None
for start_pos in range(0, longest_line_length, COLS):
  should_skip = not any(
    line[start_pos:start_pos + COLS].strip().replace(".", "")
    for line in lines[1:])
  if should_skip and skip_start is None:
    skip_start = start_pos
  elif not should_skip and skip_start is not None:
    print("Skipped %s..%s" % (
      skip_start, start_pos - 1))
    skip_start = None

  if should_skip:
    continue

  print("\n%s..%s" % (start_pos, start_pos + COLS))
  first_line = None
  for i, line in enumerate(lines):
    to_print = line[start_pos:start_pos + COLS]

    if i == 0:
      first_line = to_print
    else:
      def match(canonical, found):
        return canonical == found or found == "." or found.isdigit()
      
      to_print = "".join(
        "%s%s%s" % (
          "" if match(canonical, found) else COLOR_RED,
          found,
          "" if match(canonical, found) else COLOR_END)
        for canonical, found in zip(first_line, to_print))

    print(to_print)

if skip_start is not None:
  print("Skipped %s..%s" % (
    skip_start, longest_line_length - 1))
