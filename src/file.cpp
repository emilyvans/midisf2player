#include "file.hpp"
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stacktrace>
#include <vector>

#define panic(message)                                                         \
	{                                                                          \
		std::cerr << "PANIC: \"" << message << "\" in " << __FILE__ << ":"     \
		          << __LINE__ << "\n";                                         \
		exit(1);                                                               \
	}

File::File(std::string file_name, std::ios::openmode mode) {
	std::ifstream file(file_name, mode);

	m_cursor = 0;
	file.seekg(0, file.end);
	m_file_length = file.tellg();
	file.seekg(0, file.beg);
	m_file_bytes = new uint8_t[m_file_length];
	file.read((char *)m_file_bytes, m_file_length);
}

File::~File() { delete[] m_file_bytes; }

void File::add_to_cursor(uint32_t number) { m_cursor += number; }
uint32_t File::get_cursor() { return m_cursor; }
void File::set_cursor(uint32_t cursor) { m_cursor = cursor; }

uint8_t File::read8() {
	if (m_cursor + 1 > m_file_length)
		panic("hit end of buffer. length: " << m_file_length
		                                    << ", cursor: " << m_cursor);
	return m_file_bytes[m_cursor++];
}

uint16_t File::read16() {
	if (m_cursor + 2 > m_file_length)
		panic("hit end of buffer. length: " << m_file_length
		                                    << ", cursor: " << m_cursor);
	uint8_t byte1 = m_file_bytes[m_cursor++];
	uint8_t byte2 = m_file_bytes[m_cursor++];
	return byte2 << 8 | byte1;
}

uint32_t File::read32() {
	if (m_cursor + 4 > m_file_length)
		panic("hit end of buffer. length: " << m_file_length
		                                    << ", cursor: " << m_cursor);
	uint8_t byte1 = m_file_bytes[m_cursor++];
	uint8_t byte2 = m_file_bytes[m_cursor++];
	uint8_t byte3 = m_file_bytes[m_cursor++];
	uint8_t byte4 = m_file_bytes[m_cursor++];
	return byte4 << 24 | byte3 << 16 | byte2 << 8 | byte1;
}

uint64_t File::read64() {
	if (m_cursor + 8 > m_file_length)
		panic("hit end of buffer. length: " << m_file_length
		                                    << ", cursor: " << m_cursor);
	uint8_t byte1 = m_file_bytes[m_cursor++];
	uint8_t byte2 = m_file_bytes[m_cursor++];
	uint8_t byte3 = m_file_bytes[m_cursor++];
	uint8_t byte4 = m_file_bytes[m_cursor++];
	uint8_t byte5 = m_file_bytes[m_cursor++];
	uint8_t byte6 = m_file_bytes[m_cursor++];
	uint8_t byte7 = m_file_bytes[m_cursor++];
	uint8_t byte8 = m_file_bytes[m_cursor++];
	return (uint64_t)byte8 << 54 | (uint64_t)byte7 << 48 |
	       (uint64_t)byte6 << 40 | (uint64_t)byte5 << 32 | byte4 << 24 |
	       byte3 << 16 | byte2 << 8 | byte1;
}

std::vector<uint8_t> File::readn(uint32_t count, bool forward) {
	if (m_cursor + count > m_file_length)
		panic("hit end of buffer. length: " << m_file_length
		                                    << ", cursor: " << m_cursor);
	std::vector<uint8_t> bytes;

	if (forward) {
		for (int i = 0; i < count; i++) {
			bytes.push_back(m_file_bytes[m_cursor + i]);
		}
	} else {
		for (int i = count; i >= 0; i--) {
			bytes.push_back(m_file_bytes[m_cursor + i]);
		}
	}
	m_cursor += count;

	return bytes;
}
