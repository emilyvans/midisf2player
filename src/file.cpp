#include "file.hpp"
#include <cstdint>
#include <vector>

MyFile::MyFile(std::string file_name, std::ios::openmode mode) {}

uint8_t MyFile::read8() { return m_file_bytes[m_cursor++]; }

uint16_t MyFile::read16() {
	uint8_t byte1 = m_file_bytes[m_cursor++];
	uint8_t byte2 = m_file_bytes[m_cursor++];
	return byte2 << 8 | byte1;
}

uint32_t MyFile::read32() {
	uint8_t byte1 = m_file_bytes[m_cursor++];
	uint8_t byte2 = m_file_bytes[m_cursor++];
	uint8_t byte3 = m_file_bytes[m_cursor++];
	uint8_t byte4 = m_file_bytes[m_cursor++];
	return byte4 << 24 | byte3 << 16 | byte2 << 8 | byte1;
}

uint64_t MyFile::read64() {
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

std::vector<uint8_t> MyFile::readn(uint32_t count, bool forward) {
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
