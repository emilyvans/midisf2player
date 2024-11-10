#ifndef FILE_HPP_
#define FILE_HPP_

#include <cstdint>
#include <ios>
#include <vector>

class MyFile {

  public:
	MyFile(std::string file_name,
	       std::ios::openmode mode = std::ios::in | std::ios::binary);
	uint8_t read8();
	uint16_t read16();
	uint32_t read32();
	uint64_t read64();
	std::vector<uint8_t> readn(uint32_t count, bool forward = false);

  private:
  public:
  private:
	uint32_t m_cursor;
	uint32_t m_file_length;
	uint8_t *m_file_bytes;
};

#endif // FILE_HPP_
