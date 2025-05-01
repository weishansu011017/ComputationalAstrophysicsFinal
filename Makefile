#=== Toolchain ===#
CXX      ?= clang++

#=== Path flags ===#
INCLUDE_DIRS := -Iinclude -Isrc
LIB_DIRS     := -L/opt/homebrew/lib

#=== Feature flags ===#


#=== External library flags ===#
EXTERNAL_LIBS := hdf5

#=== Aggregate final flags ===#
CXXFLAGS  := -std=c++17 $(INCLUDE_DIRS) $(shell pkg-config --cflags $(EXTERNAL_LIBS))
LDFLAGS := $(LIB_DIRS) $(shell pkg-config --libs $(EXTERNAL_LIBS))

#=== File structure ===#
SRC_DIR     := src
TEST_DIR    := test
BUILD_DIR   := build

#=== Searching all test/*.cpp ===#
TESTS := $(basename $(notdir $(wildcard $(TEST_DIR)/*.cpp)))

#=== Target of testfiles ===#
TESTSTARGETS := $(addprefix $(BUILD_DIR)/test/, $(TESTS))

#=== Default make ===#
all:
	@echo "### Compiler        : $(CXX)"
	@echo "### External libs   : $(EXTERNAL_LIBS)"
	@echo "Tip: If your dependencies are not installed in system paths,"
	@echo "        you may need to set:"
	@echo "        export PKG_CONFIG_PATH=/your/custom/path/lib/pkgconfig"
	@echo "Default build: main simulation target is not yet implemented."

#=== Compiling TestXXXXX.cpp to build/test/ ===#
##=== General rule for specific test code folder to test/
$(BUILD_DIR)/test/%: $(TEST_DIR)/%.cpp $(SRC_DIR)/*.cpp
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