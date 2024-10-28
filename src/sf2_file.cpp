#include "sf2_file.hpp"
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

uint8_t read8(uint8_t *bytes, uint32_t file_length, uint32_t &cursor);
uint16_t read16(uint8_t *bytes, uint32_t file_length, uint32_t &cursor);
uint32_t read32(uint8_t *bytes, uint32_t file_length, uint32_t &cursor);
std::vector<uint8_t> read_number(uint8_t *in_bytes, uint32_t file_length,
                                 uint32_t &cursor, uint32_t count);
std::vector<uint8_t> read_number_forward(uint8_t *in_bytes,
                                         uint32_t file_length, uint32_t &cursor,
                                         uint32_t count);

SF2_FILE::SF2_FILE(std::string file_name) {
	std::ifstream file(file_name, std::ios::binary);
	file.seekg(0, file.end);
	file_length = file.tellg();
	file.seekg(0, file.beg);
	uint8_t *file_bytes = new uint8_t[file_length];
	uint32_t cursor = 0;
	file.read((char *)file_bytes, file_length);

	std::function<uint8_t()> read1 = [file_bytes, this, &cursor]() -> uint8_t {
		return read8(file_bytes, file_length, cursor);
	};
	std::function<uint16_t()> read2 = [file_bytes, this,
	                                   &cursor]() -> uint16_t {
		return read16(file_bytes, file_length, cursor);
	};
	std::function<uint32_t()> read4 = [file_bytes, this,
	                                   &cursor]() -> uint32_t {
		return read32(file_bytes, file_length, cursor);
	};
	std::function<std::vector<uint8_t>(uint32_t count)> readn =
	    [file_bytes, this, &cursor](uint32_t count) {
		    return read_number(file_bytes, file_length, cursor, count);
	    };
	std::function<std::vector<uint8_t>(uint32_t count)> readn_forward =
	    [file_bytes, this, &cursor](uint32_t count) {
		    return read_number_forward(file_bytes, file_length, cursor, count);
	    };

	read4(); // RIFF
	read4(); // riff chunk length
	read4(); // SFBK

	// INFO
	read4(); // LIST
	read4(); // LIST chunk length
	read4(); // INFO
	read4(); // ifil
	read4(); // ifil length
	major_version = read2();
	minor_version = read2();
	read4();                        // isng
	uint32_t isng_length = read4(); // isng length
	sound_engine.assign((char *)readn_forward(isng_length).data(), isng_length);
	read4(); // INAM
	uint32_t inam_length = read4();
	name.assign((char *)readn_forward(inam_length).data(), inam_length);
	read4(); // ICRD
	uint32_t icrd_length = read4();
	date_of_creation.assign((char *)readn_forward(icrd_length).data(),
	                        icrd_length);
	read4(); // ISFT
	uint32_t isft_length = read4();
	tool.assign((char *)readn_forward(isft_length).data(), isft_length);

	// sdta
	read4(); // LIST
	read4(); // LIST chunk length
	read4(); // sdta
	read4(); // smpl
	uint32_t smpl_length = read4();

	samples = readn_forward(smpl_length);
	std::ofstream audio("audio.bin", std::ios::binary);

	for (uint8_t smpl : samples) {
		audio << smpl;
	}

	audio.close();

	// pdta
	read4();                        // LIST
	read4();                        // LIST chunk length
	read4();                        // pdta
	read4();                        // phdr
	uint32_t phdr_length = read4(); // phdr chunk length

	for (int i = 0; i < phdr_length; i += 38) {
		sf_preset_header preset_header;
		std::vector<uint8_t> preset_name = readn_forward(20);
		memcpy(preset_header.preset_name, preset_name.data(), 20);
		preset_header.preset = read2();
		preset_header.bank = read2();
		preset_header.perset_bag_ndx = read2();
		preset_header.library = read4();
		preset_header.genre = read4();
		preset_header.morphology = read4();
		preset_headers.push_back(preset_header);
	}

	std::cout << std::hex << read4() << "\n";

	delete[] file_bytes;
}

uint8_t read8(uint8_t *bytes, uint32_t file_length, uint32_t &cursor) {
	if (file_length < cursor + 1) {
		std::cerr << "cursor(" << cursor + 1 << ") is bigger than file_length("
		          << file_length << ")\n";
		assert(false);
	}
	return bytes[cursor++];
}

uint16_t read16(uint8_t *bytes, uint32_t file_length, uint32_t &cursor) {
	if (file_length < cursor + 2) {
		std::cerr << "cursor(" << cursor + 2 << ") is bigger than file_length("
		          << file_length << ")\n";
		assert(false);
	}
	uint8_t byte1 = bytes[cursor++];
	uint8_t byte2 = bytes[cursor++];

	return byte2 << 8 | byte1;
}

uint32_t read32(uint8_t *bytes, uint32_t file_length, uint32_t &cursor) {
	if (file_length < cursor + 4) {
		std::cerr << "cursor(" << cursor + 4 << ") is bigger than file_length("
		          << file_length << ")\n";
		assert(false);
	}
	uint8_t byte1 = bytes[cursor++];
	uint8_t byte2 = bytes[cursor++];
	uint8_t byte3 = bytes[cursor++];
	uint8_t byte4 = bytes[cursor++];

	return byte4 << 24 | byte3 << 16 | byte2 << 8 | byte1;
}

std::vector<uint8_t> read_number(uint8_t *in_bytes, uint32_t file_length,
                                 uint32_t &cursor, uint32_t count) {
	if (file_length < cursor + count) {
		std::cerr << "cursor(" << cursor + count
		          << ") is bigger than file_length(" << file_length << ")\n";
		assert(false);
	}
	int index = count - 1;
	uint8_t *out_bytes = new uint8_t[count];
	while (index >= 0) {
		out_bytes[index] = in_bytes[cursor++];
		index -= 1;
	}

	std::vector<uint8_t> out(out_bytes, out_bytes + count);
	return out;
}

std::vector<uint8_t> read_number_forward(uint8_t *in_bytes,
                                         uint32_t file_length, uint32_t &cursor,
                                         uint32_t count) {
	if (file_length < cursor + count) {
		std::cerr << "cursor(" << cursor + count
		          << ") is bigger than file_length(" << file_length << ")\n";
		assert(false);
	}
	std::vector<uint8_t> out;

	for (int i = 0; i < count; i++) {
		out.push_back(in_bytes[cursor++]);
	}

	return out;
}
