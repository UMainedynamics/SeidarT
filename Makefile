FC = gfortran
FFLAGS = -J./mod_dir -O2 -fPIC
SRC = src/seidart/fortran/cpmlfdtd.f95
OBJ_DIR = build
MOD_DIR = mod_dir
OBJ = $(OBJ_DIR)/cpmlfdtd.o
MOD = $(MOD_DIR)/cpmlfdtd.mod
EXT = $(OBJ_DIR)/cpmlfdtd.so

.PHONY: all clean  # Force these targets to always run

all: $(EXT)
	@echo "Build complete!"

# Ensure the build directories exist
$(OBJ_DIR) $(MOD_DIR):
	@echo "Creating directories..."
	mkdir -p $(OBJ_DIR) $(MOD_DIR)

# Compile the Fortran source and generate .mod and .o files
$(OBJ): $(SRC) $(OBJ_DIR) $(MOD_DIR)
	@echo "Compiling Fortran source..."
	$(FC) $(FFLAGS) -c $(SRC) -o $(OBJ)

# Link the Fortran object file into a shared library for Python
$(EXT): $(OBJ)
	@echo "Linking object file into shared library..."
	$(FC) -shared -o $(EXT) $(OBJ)

clean:
	@echo "Cleaning up..."
	rm -rf $(OBJ_DIR) $(MOD_DIR)
