#include <fstream>
#include <iostream>
#include <optional>
#include <vector>
#include <bitset>

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

	return byte4 << 24 | byte3 << 16 | byte2 << 8 | byte1;
}

uint8_t *readn_little_endian(std::ifstream &file, int n) {
	char *result = new char[](n);
	file.read(result, n);

	uint8_t *result_swap = new uint8_t[](n);

	for (int i = sizeof(result), i_swap = 0; i > 0; i--, i_swap++) {
		result_swap[i_swap] = result[i];
	}

	return result_swap;
}

struct WAVE_FILE {
	char magic_bytes_riff[4]; // "riff"
	char magic_bytes[4];      // "WAVE"
	uint32_t file_length;
	uint16_t format_tag;        // format category
	uint16_t channels;          // number of channels
	uint32_t samples_per_sec;   // sampling rate
	uint32_t avg_bytes_per_sec; // for buffer estimation
	uint16_t block_align;       // data block size
	uint16_t bits_per_sample;   // the bits per sample
	uint32_t data_length;
	std::vector<double> samples;
};

std::optional<WAVE_FILE> read_wave_file(const char *file_path) {
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
		std::cout << "format isn't 1 but " << wave.format_tag << "\n";
		return std::nullopt;
	}
	wave.channels = read2(file);
	wave.samples_per_sec = read4(file);
	wave.avg_bytes_per_sec = read4(file);
	wave.block_align = read2(file);
	wave.bits_per_sample = read2(file);

	read4(file); // "data"
	wave.data_length = read4(file);
	int sample_increment = ceilf(wave.bits_per_sample / (float)8);
	int64_t max = 255;
	int64_t min = 0;
	if (wave.bits_per_sample > 8) {
		int range = pow(2, wave.bits_per_sample);
		max = range / 2 - 1;
		min = -1 * (range / 2);
	}
	std::cout << "min: " << min << " max: " << max << "\n";
	for (int i = 0; i < sample_increment * 2; i += sample_increment) {
		uint8_t *raw_sample_bytes = readn_little_endian(file, sample_increment);
		uint64_t raw_sample = 0;
		int64_t sample_int = 0;
		float sample = 0;
		for (int i = 0; i < sample_increment; i++) {
			raw_sample = raw_sample << 8 | raw_sample_bytes[i];
		}
		bool isNegative = false;
		
		std::cout << std::bitset<16>(raw_sample) << " ";

		if ((raw_sample & ((uint64_t)0b1 << sample_increment * 8 - 1)) >= 1) {
			isNegative = true;
		}

		std::cout << " " << (isNegative ? "true" : "false") << "\n"; 
		raw_sample >>= (sample_increment * 8) - wave.bits_per_sample;

		if (isNegative) {
			sample_int = ~raw_sample - min;
		} else {
			sample_int = raw_sample;
		}

		sample = (float)sample_int / (float)max;

		std::cout << "bytes: " << std::hex << (uint16_t *)raw_sample
							<< " raw: " << (uint64_t)raw_sample << " : " << std::dec
							<< (uint64_t)raw_sample << " : act: " << sample
							<< " bin: " << std::bitset<16>(raw_sample) << " : "
							<< wave.bits_per_sample << " : " << isNegative 
							<< "\n";
	}

	return std::make_optional(wave);
}

int main() {
	std::optional<WAVE_FILE> wave_opt = read_wave_file("test.wav");
	if (!wave_opt.has_value()) {
		std::cout
				<< "coudn't parse wave file. it either doesn't exist or isn't format 1";
		std::cin.get();
		return -1;
	}
	WAVE_FILE wave = wave_opt.value();

	std::cin.get();
}