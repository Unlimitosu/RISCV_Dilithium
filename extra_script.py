Import("env")

# General options that are passed to the C and C++ compilers
env.Append(CCFLAGS=["-O3", ""])