#include <iostream>
#include <fstream>
#include <optional>
#include <vector>

uint8_t read1(std::ifstream &file) { 
	char result[1];

	file.read(result, sizeof(result));

	return result[0];
}

uint16_t read2(std::ifstream &file) {
	char result[2];
	file.read(result, sizeof(result));

	uint8_t byte1 = result[0];
	uint8_t byte2 = result[1];

	return byte2 << 8 | byte1;
}

uint32_t read3(std::ifstream &file) {
	char result[3];
	file.read(result, sizeof(result));

	uint8_t byte1 = result[0];
	uint8_t byte2 = result[1];
	uint8_t byte3 = result[2];

	return byte3 << 16 | byte2 << 8 | byte3;
}

uint32_t read4(std::ifstream &file) {
	char result[4];
	file.read(result, sizeof(result));

	uint8_t byte1 = result[0];
	uint8_t byte2 = result[1];
	uint8_t byte3 = result[2];
	uint8_t byte4 = result[3];

	return byte1 << 24 | byte2 << 16 | byte3 << 8 | byte4;
}

union sound_sample {

};

struct WAVE_FILE{
	char magic_bytes_riff[4];	// "riff"
	char magic_bytes[4];		// "WAVE"
	uint32_t file_length;
	uint16_t format_tag;		// format category
	uint16_t channels;			// number of channels
	uint32_t samples_per_sec;	// sampling rate
	uint32_t avg_bytes_per_sec;	// for buffer estimation
	uint16_t block_align;		// data block size
	uint16_t bits_per_sample;	// the bits per sample
	uint32_t data_length;
	union {
		std::vector<int8_t> int8data;
		std::vector<int16_t> int16data;
		std::vector<int32_t> int32data;
		std::vector<uint8_t> uint8data;
		std::vector<uint16_t> uint16data;
		std::vector<uint32_t> uint32data;
		std::vector<uint8_t> bytes;
	};
};

std::optional<WAVE_FILE> read_wave_file(const char* file_path) {
	std::ifstream file(file_path, std::ios::binary);
	if (!file.good()) {
		return std::nullopt;
	}
	WAVE_FILE wave;
	file.read(wave.magic_bytes_riff, 4);
	wave.file_length = read4(file);
	file.read(wave.magic_bytes, 4);
	read4(file); // "fmt "
	read4(file); // fmt chunk sizes
	wave.format_tag = read2(file);
	if (wave.format_tag != 1) {
		std::cout << "format isn't 1 but " << wave.format_tag
							<< "\n";
		return std::nullopt;
	}
	wave.channels = read2(file);
	wave.samples_per_sec = read4(file);
	wave.avg_bytes_per_sec = read4(file);
	wave.block_align = read2(file);
	wave.bits_per_sample = read2(file);

	read4(file); // "data"
	wave.data_length = read4(file);



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