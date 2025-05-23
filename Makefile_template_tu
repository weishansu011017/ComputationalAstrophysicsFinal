#=== Toolchain ===#
LLVM_DIR := /opt/homebrew/opt/llvm
CXX := g++

#=== Path flags ===#
INCLUDE_DIRS := -Iinclude -Isrc -Iutil/tomlplusplus/include
LIB_DIRS     := -L/usr/lib/gcc/x86_64-linux-gnu/11/

#=== Feature flags ===#


#=== External library flags ===#
EXTERNAL_LIBS := hdf5

#=== Aggregate final flags ===#
CXXFLAGS := -std=c++17 -O3 -march=native -fopenmp $(INCLUDE_DIRS) \
            $(shell pkg-config --cflags $(EXTERNAL_LIBS))
LDFLAGS := $(LIB_DIRS) $(shell pkg-config --libs $(EXTERNAL_LIBS)) \
           -L$(LLVM_DIR)/lib -Wl,-rpath,$(LLVM_DIR)/lib

#=== File structure ===#
CORE_DIR     := src/core
MAIN_DIR     := src/main
TEST_DIR    := test
BUILD_DIR   := build

#=== Searching all test/*.cpp ===#
TESTS := $(basename $(notdir $(wildcard $(TEST_DIR)/*.cpp)))

#=== Target of testfiles ===#
TESTSTARGETS := $(addprefix $(BUILD_DIR)/test/, $(TESTS))

#=== Target of setup.cpp ===#
SETUPTARGET := setup
SETUP_SRC   := $(MAIN_DIR)/setup.cpp $(wildcard $(CORE_DIR)/*.cpp)

#=== Default make ===#
all:
	@echo "### Compiler        : $(CXX)"
	@echo "### External libs   : $(EXTERNAL_LIBS)"
	@echo "Tip: If your dependencies are not installed in system paths,"
	@echo "        you may need to set:"
	@echo "        export PKG_CONFIG_PATH=/your/custom/path/lib/pkgconfig"
	@echo "Default build: main simulation target is not yet implemented."

#=== Default Initial Condition constructor ===#
setup:
$(SETUPTARGET): $(SETUP_SRC)
	@echo "### Compiler        : $(CXX)"
	@echo "### External libs   : $(EXTERNAL_LIBS)"
	@echo "Tip: If your dependencies are not installed in system paths,"
	@echo "        you may need to set:"
	@echo "        export PKG_CONFIG_PATH=/your/custom/path/lib/pkgconfig"
	@echo "                   "
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

#=== Compiling TestXXXXX.cpp to build/test/ ===#
##=== General rule for specific test code folder to test/
$(BUILD_DIR)/test/%: $(TEST_DIR)/%.cpp $(CORE_DIR)/*.cpp
	@mkdir -p $(BUILD_DIR)/test
	$(CXX) $(CXXFLAGS) $(filter %.cpp, $^) -o $@ $(LDFLAGS)
##=== Compiling specific TestXXXXX.cpp ===##
$(TESTS): %: $(BUILD_DIR)/test/%
	@echo "Built $@"
##=== Compiling all TestXXXXX.cpp to build/ ===##
testall:
	@echo "[Test Build] Building $(words $(TESTSTARGETS)) test programs..."
	@$(MAKE) --no-print-directory $(TESTSTARGETS)
	@echo "All tests built successfully!"

#=== Clean all build ===#
clean:
	@find $(BUILD_DIR) -type f ! -name ".gitkeep" -delete
	@find $(BUILD_DIR) -type d -empty -not -path "$(BUILD_DIR)" -delete
	@echo "Cleaned build/"