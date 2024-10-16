cmake_minimum_required(VERSION 3.29.3)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 23)

project(fft_test)
project(audio_test)
project(wav_test)

add_executable(fft_test src/main.cpp)
add_executable(audio_test src/main_from_example.cpp)
add_executable(wav_test src/main_wav.cpp)

set(BUILD_SHARED_LIBS OFF)

add_subdirectory(vendor/portaudio)

target_link_libraries(fft_test PortAudio)
target_link_libraries(audio_test PortAudio)
target_link_libraries(wav_test PortAudio)
# target_link_libraries(midi ${OPENGL_opengl_LIBRARY} glad glfw)

#file(COPY src/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
#add_custom_command(TARGET midi POST_BUILD 
#		COMMAND ${CMAKE_COMMAND} -E copy_directory
#				${CMAKE_CURRENT_SOURCE_DIR}/sample
#				${CMAKE_CURRENT_BINARY_DIR}/sample)
