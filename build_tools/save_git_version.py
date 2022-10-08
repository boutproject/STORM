#!/usr/bin/env python3

import argparse
from boututils.run_wrapper import shell_safe
from pathlib import Path
from string import ascii_uppercase

try:
    from git import Repo, InvalidGitRepositoryError
except ImportError:
    # print a return code to stdout so we can capture and test it in the makefile
    print(11)
    raise ImportError(
        "gitpython not installed. You can install with 'conda install gitpython' or "
        "'pip3 install --user gitpython'"
    )

scriptdir = Path(__file__).resolve().parent

parser = argparse.ArgumentParser()
parser.add_argument(
    "output_directory",
    type=str,
    nargs="?",
    default=scriptdir.parent.joinpath("storm3d"),
)
parser.add_argument("cmake_directory", type=str, nargs="?", default=None)
args = parser.parse_args()
output_directory = Path(args.output_directory)
cmake_directory = None if args.cmake_directory is None else Path(args.cmake_directory)

repo = Repo(scriptdir.parent)

diff = repo.git.diff()

# Get version info for storm-configs, if it is being used
try:
    storm_configs_repo = Repo(scriptdir.resolve().parent.parent)
except InvalidGitRepositoryError:
    storm_configs_repo = None

if storm_configs_repo is not None:
    # Check that storm_configs_repo is really storm-configs
    *_, first_commit = storm_configs_repo.iter_commits()
    if str(first_commit) != "f45792f72fb401014459987d1064113fa572c625":
        storm_configs_repo = None
    else:
        storm_configs_diff = storm_configs_repo.git.diff()
        bout_configs_repo = Repo(
            scriptdir.resolve().parent.parent.joinpath("BOUT-configs")
        )
        bout_configs_diff = bout_configs_repo.git.diff()

# Get contents of CMakeCache.txt, if CMake is being used
if cmake_directory is not None:
    with open(cmake_directory.joinpath("CMakeCache.txt"), "r") as f:
        storm_cmake_cache = f.read()

# Get modules, if the system has a 'module' command
try:
    _, module_list = shell_safe("module list", pipe=True)
except RuntimeError:
    module_list = None

# Find a string to use for the delimiter
if storm_configs_repo is not None:
    combined_diff = diff + "\n" + storm_configs_diff + "\n" + bout_configs_diff
else:
    combined_diff = diff
if cmake_directory is not None:
    combined_diff += "\n" + storm_cmake_cache
if module_list is not None:
    combined_diff += "\n" + module_list
delimiter = "_STORM_DIFF_"
for letter in ascii_uppercase:
    if ")" + delimiter not in combined_diff:
        break
    delimiter = "_STORM_DIFF_" + letter + "_"
if ")" + delimiter in combined_diff:
    # print a return code to stdout so we can capture and test it in the makefile
    print(12)
    raise ValueError(
        "save_git_version.py failed to find a delimiter that is not in the git diff"
    )

with open(output_directory.joinpath("storm_version.hxx"), "w") as f:
    f.write('constexpr auto storm_git_hash = "')
    f.write(repo.git.describe(abbrev=40, dirty=True, always=True, tags=True, long=True))
    f.write('";\n')
    f.write('constexpr auto storm_git_diff = R"' + delimiter + "(")
    f.write(diff)
    f.write(")" + delimiter + '";\n')

    if storm_configs_repo is None:
        f.write('constexpr auto storm_configs_git_hash = "";\n')
        f.write('constexpr auto storm_configs_git_diff = "";\n')
        f.write('constexpr auto bout_configs_git_diff = "";\n')
    else:
        f.write('constexpr auto storm_configs_git_hash = "')
        f.write(
            storm_configs_repo.git.describe(
                abbrev=40,
                dirty=True,
                always=True,
                tags=True,
            ),
        )
        f.write('";\n')
        f.write('constexpr auto storm_configs_git_diff = R"' + delimiter + "(")
        f.write(storm_configs_diff)
        f.write(")" + delimiter + '";\n')
        f.write('constexpr auto bout_configs_git_diff = R"' + delimiter + "(")
        f.write(bout_configs_diff)
        f.write(")" + delimiter + '";\n')

    if cmake_directory is None:
        f.write('std::string storm_cmake_cache = "";\n')
    else:
        f.write('std::string storm_cmake_cache = R"' + delimiter + "(")
        f.write(storm_cmake_cache)
        f.write(")" + delimiter + '";\n')

    if module_list is None:
        f.write('std::string module_list = "";\n')
    else:
        f.write('std::string module_list = R"' + delimiter + "(")
        f.write(module_list)
        f.write(")" + delimiter + '";\n')

# print a return code to stdout so we can capture and test it in the makefile
print(0)
