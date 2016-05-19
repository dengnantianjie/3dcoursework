=== UNIX ===

cd extern/glfw-3.1.2
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_SHARED_LIBS=ON -DGLFW_BUILD_DOCS=OFF -DGLFW_BUILD_EXAMPLES=OFF -DGLFW_BUILD_TESTS=OFF  -L ..
make -j 8
make install
cd ../../..
mkdir build
cd build
cmake ..
make -j 4
./meshviewer bunny000.ply

=== WINDOWS ===
- download & install cmake for windows from the cmake website
- if the installer complains that the path string is too long, add the bin directory of cmake to the path manually (for example, add 'C:\Program Files (x86)\CMake\bin' to the path)
- test that cmake works by opening a command prompt and running 'cmake', cmake should display the usage info.
- extract the framework and go to the 'lib' directory of the framework
- execute the follwing commands:
cd extern/glfw-3.1.2
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_SHARED_LIBS=ON -DGLFW_BUILD_DOCS=OFF -DGLFW_BUILD_EXAMPLES=OFF -DGLFW_BUILD_TESTS=OFF -DGLFW_DLL=ON -G "Visual Studio 12 2013 Win64" ..
(or any other supported target instead of "Visual Studio 12 2013", type in an unsupported target to see a list of all targets)
cmake --build . --config Release --target install
cd ../../..
mkdir build
cd build
cmake -G "Visual Studio 12 2013 Win64" ..
cmake --build . --config Release
copy ..\extern\glew-1.13.0\bin\Release\x64\glew32.dll .\Release\glew32.dll
copy ..\extern\glfw-3.1.2\install\lib\glfw3.dll .\Release\glfw3.dll
cd Release
meshviewer bun000.ply
