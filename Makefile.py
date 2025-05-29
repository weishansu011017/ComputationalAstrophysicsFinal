### Generate by Chat GPT-o3
#!/usr/bin/env python
# build.py  --- Windows 取代 Makefile 的簡易 build script
#
# 用法：
#   python build.py                # same as “make simulation”
#   python build.py simulation     # build simulation.exe
#   python build.py setup          # build setup.exe
#   python build.py testall        # build every test/*.cpp
#   python build.py clean          # remove build artefacts
#
# 如要修改 HDF5 位置、編譯器或其他 flag，直接編輯下方常數即可。

import sys
import subprocess
import os
from pathlib import Path
from shutil import rmtree
from glob import glob
import shlex

# ──────────────────────────────────────────────────────────────
# 🔧 1. 路徑與工具設定
# ──────────────────────────────────────────────────────────────
# 若你用的是 LLVM for Windows（例如 scoop / choco 裝的），通常只要把 clang++
# 加進 PATH 即可；如果想強迫使用 LLVM 安裝資料夾，取消註解 LLVM_DIR 再下方用它。
# LLVM_DIR = Path(r"C:\Program Files\LLVM")        # 例：LLVM 安裝根目錄
CXX = r"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\Llvm\x64\bin\clang-cl.exe"                                     # or "g++", "cl"

# HDF5 預設安裝路徑（請依自己的版本調整）
HDF5_ROOT = Path(r"C:\Program Files\HDF_Group\HDF5\1.14.6")
HDF5_INC  = HDF5_ROOT / "include"
HDF5_LIB  = HDF5_ROOT / "lib"
SZIP_ROOT = Path("E:/programming/libaec/build/src/Release")

# 專案結構
CORE_DIR  = Path("src/core")
MAIN_DIR  = Path("src/main")
TEST_DIR  = Path("test")
BUILD_DIR = Path("build")
(BUILD_DIR / "test").mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────
# 🔧 2. 編譯／連結旗標
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

# ⚠️ Windows 不支援 -march=native 給 MSVC；若你改用 cl.exe 請刪掉。
COMMON_CXXFLAGS = [
    "/std:c++17",
    "/O2",
    "/MD",        # ✅ 配合 MSVC 預設 CRT（HDF5 預設用 /MD 編譯）
    "/openmp",
    "/EHsc",      # exception handler
    "/bigobj",
    "-DM_PI=3.14159265358979323846",
]


# 直接列絕對路徑，避免空白路徑被 linker 誤判
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
# 🔧 3. 公用函式
# ──────────────────────────────────────────────────────────────
def run(cmd: list[str]) -> None:
    """執行編譯指令並在失敗時終止腳本。"""
    print(">>", " ".join(shlex.quote(a) for a in cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)

def compile_exe(output: Path, sources: list[Path]) -> None:
    """將 sources 編譯成指定 exe。"""
    cmd = [CXX, *COMMON_CXXFLAGS, *INCLUDE_DIRS,
           *map(str, sources),
           "-o", str(output),
           *LIB_DIRS, *LDFLAGS]
    run(cmd)

# ──────────────────────────────────────────────────────────────
# 🔧 4. 目標邏輯
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
# 🔧 5. 主程式
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
