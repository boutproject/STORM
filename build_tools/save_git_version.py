#!/usr/bin/env python3

from pathlib import Path
from string import ascii_uppercase

try:
    from git import Repo
except ImportError:
    # print a return code to stdout so we can capture and test it in the makefile
    print(11)
    raise

scriptdir = Path(__file__).parent

repo = Repo(scriptdir.parent)

diff = repo.git.diff()

# Find a string to use for the delimiter
delimiter = "_STORM_DIFF_"
for letter in ascii_uppercase:
    if ")" + delimiter not in diff:
        delimiter_success = True
        break
    delimiter = "_STORM_DIFF_" + letter + "_"
if ")" + delimiter in diff:
    # print a return code to stdout so we can capture and test it in the makefile
    print(12)
    raise ValueError(
        "save_git_version.py failed to find a delimiter that is not in the git diff"
    )

with open("storm_version.hxx", "w") as f:
    f.write("constexpr auto storm_git_hash = \"")
    f.write(repo.git.describe(abbrev=40, dirty=True, always=True, tags=True))
    f.write("\";\n")
    f.write("constexpr auto storm_git_diff = R\"" + delimiter + "(")
    f.write(diff)
    f.write(")" + delimiter + "\";\n")

# print a return code to stdout so we can capture and test it in the makefile
print(0)
