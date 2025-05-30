### Generate by Chat GPT-o3
#!/usr/bin/env python
# build.py  --- Windows 取代 Makefile 的簡易 build script
#
# Usage:
#   python build.py                # same as “make simulation”
#   python build.py simulation     # build simulation.exe
#   python build.py setup          # build setup.exe
#   python build.py testall        # build every test/*.cpp
#   python build.py clean          # remove build artefacts
#

import sys
import subprocess
import os
from pathlib import Path
from shutil import rmtree
from glob import glob
import shlex

# ──────────────────────────────────────────────────────────────
#       Path & Tool chain
# ──────────────────────────────────────────────────────────────
CXX = r"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\Llvm\x64\bin\clang-cl.exe"                        

# HDF5 Path
HDF5_ROOT = Path(r"C:\Program Files\HDF_Group\HDF5\1.14.6")
HDF5_INC  = HDF5_ROOT / "include"
HDF5_LIB  = HDF5_ROOT / "lib"
SZIP_ROOT = Path("E:/programming/libaec/build/src/Release")

# Project Path
CORE_DIR  = Path("src/core")
MAIN_DIR  = Path("src/main")
TEST_DIR  = Path("test")
BUILD_DIR = Path("build")
(BUILD_DIR / "test").mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────
#       Flag of compiling
# ──────────────────────────────────────────────────────────────
INCLUDE_DIRS = [
    "/Iinclude",
    "/Isrc",
    "/Iutil/tomlplusplus/include",
    f"/I{HDF5_INC}"
]

LIB_DIRS = [
    f"/link",
    "/machine:x64",
    f"/LIBPATH:{HDF5_LIB}"
]

COMMON_CXXFLAGS = [
    "/std:c++17",
    "/Ox",
    "/MD",        
    "/openmp",
    "/EHsc",    
    "/bigobj",
    "-DM_PI=3.14159265358979323846",
]


LIBS = [
    HDF5_LIB / "libhdf5_cpp.lib",
    HDF5_LIB / "libhdf5.lib",
    HDF5_LIB / "libhdf5_hl.lib",
    HDF5_LIB / "zlib-static.lib",     
    SZIP_ROOT / "aec-static.lib",
    SZIP_ROOT / "szip-static.lib",
    "shlwapi.lib",
]
LDFLAGS = [str(lib) for lib in LIBS]

# ──────────────────────────────────────────────────────────────
#       Common Method
# ──────────────────────────────────────────────────────────────
def run(cmd: list[str]) -> None:
    """Execute the compiling command"""
    print(">>", " ".join(shlex.quote(a) for a in cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)

def compile_exe(output: Path, sources: list[Path]) -> None:
    """Compile program into .exe"""
    cmd = [CXX, *COMMON_CXXFLAGS, *INCLUDE_DIRS,
           *map(str, sources),
           "-o", str(output),
           *LIB_DIRS, *LDFLAGS]
    run(cmd)

# ──────────────────────────────────────────────────────────────
#           Target
# ──────────────────────────────────────────────────────────────
def target_simulation() -> None:
    core = list(Path(CORE_DIR).glob("*.cpp"))
    src  = [MAIN_DIR / "simulation.cpp", *core]
    compile_exe("simulation.exe", src)

def target_setup() -> None:
    core = list(Path(CORE_DIR).glob("*.cpp"))
    src  = [MAIN_DIR / "setup.cpp", *core]
    compile_exe("setup.exe", src)

def target_testall() -> None:
    tests = sorted(Path(TEST_DIR).glob("*.cpp"))
    if not tests:
        print("No test/*.cpp found.")
        return
    print(f"[Test Build] Building {len(tests)} test programs …")
    core = list(Path(CORE_DIR).glob("*.cpp"))
    for tcpp in tests:
        exe = BUILD_DIR / "test" / f"{tcpp.stem}.exe"
        compile_exe(exe, [tcpp, *core])
    print("All tests built successfully!")

def target_clean() -> None:
    if BUILD_DIR.exists():
        for item in BUILD_DIR.iterdir():
            if item.name == ".gitkeep":
                continue
            if item.is_dir():
                rmtree(item)
            else:
                item.unlink()
        print("Cleaned build/")
    else:
        print("Nothing to clean.")

# ──────────────────────────────────────────────────────────────
#       Main
# ──────────────────────────────────────────────────────────────
DEFAULT_TARGET = "simulation"
TARGETS = {
    "simulation": target_simulation,
    "setup"     : target_setup,
    "testall"   : target_testall,
    "clean"     : target_clean,
}

if __name__ == "__main__":
    tgt = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_TARGET
    if tgt not in TARGETS:
        print("Unknown target. Use one of:", ", ".join(TARGETS))
        sys.exit(1)
    TARGETS[tgt]()
