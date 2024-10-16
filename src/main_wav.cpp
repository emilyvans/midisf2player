#include <iostream>
#include <fstream>
#include <optional>

uint8_t read1(std::ifstream &file) { 
	char result[1];

	file.read(result, sizeof(result));

	return result[0];
}

uint16_t read2(std::ifstream &file) {
	char result[2]; // 0x01 0x00

	file.read(result, sizeof(result));

	return result[0] << 8 | result[1];
}

uint32_t read4(std::ifstream &file) {
	char result[4];

	file.read(result, sizeof(result));

	return result[0] << 24 | result[0] << 16 | result[0] << 8 | result[3];
}

struct WAVE_FILE{
	char magic_bytes[4];		// "WAVE"
	uint16_t format_tag;		// format category
	uint16_t channels;			// number of channels
	uint32_t samples_per_sec;	// sampling rate
	uint32_t avg_bytes_per_sec;	// for buffer estimation
	uint16_t block_align;		// data block size
	uint16_t bits_per_sample;	// the bits per sample
};

std::optional<WAVE_FILE> read_wave_file(const char* file_path) {
	std::ifstream file("test.wav", std::ios::binary);
	if (!file.good()) {
		return std::nullopt;
	}
	WAVE_FILE wave;
	file.read(wave.magic_bytes, 4);
	wave.format_tag = read1(file);
	if (wave.format_tag != 1) {
		return std::nullopt;
	}
	wave.channels = read1(file);

	return std::make_optional(wave);
}

int main() { 
	std::optional<WAVE_FILE> wave_opt = read_wave_file("test.wav");
	if (!wave_opt.has_value()) {
		std::cout << "coudn't parse wave file. it either doesn't exist or isn't format 1";
		std::cin.get();
		return -1;
	}
	WAVE_FILE wave = wave_opt
	                     .value();

	std::cin.get();
}