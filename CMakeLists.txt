cmake_minimum_required(VERSION 3.29.3)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(portAudioTest)

add_executable(portAudioTest src/main.cpp)

set(BUILD_SHARED_LIBS ON)

add_subdirectory(vendor/portaudio)

target_link_libraries(portAudioTest portaudio)
# target_link_libraries(midi ${OPENGL_opengl_LIBRARY} glad glfw)

#file(COPY src/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
#add_custom_command(TARGET midi POST_BUILD 
#		COMMAND ${CMAKE_COMMAND} -E copy_directory
#				${CMAKE_CURRENT_SOURCE_DIR}/sample
#				${CMAKE_CURRENT_BINARY_DIR}/sample)
