#include "sf2_file.hpp"
#include "file.hpp"
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

SF2File::SF2File(std::string file_name) {
	File file(file_name);

	file.read32(); // RIFF
	file.read32(); // riff chunk length
	file.read32(); // SFBK

	// INFO
	file.read32(); // LIST
	file.read32(); // LIST chunk length
	file.read32(); // INFO
	file.read32(); // ifil
	file.read32(); // ifil length
	major_version = file.read16();
	minor_version = file.read16();
	file.read32();                        // isng
	uint32_t isng_length = file.read32(); // isng length
	sound_engine.assign((char *)file.readn(isng_length, true).data(),
	                    isng_length);
	file.read32(); // INAM
	uint32_t inam_length = file.read32();
	name.assign((char *)file.readn(inam_length, true).data(), inam_length);
	file.read32(); // ICRD
	uint32_t icrd_length = file.read32();
	date_of_creation.assign((char *)file.readn(icrd_length, true).data(),
	                        icrd_length);
	file.read32(); // ISFT
	uint32_t isft_length = file.read32();
	tool.assign((char *)file.readn(isft_length, true).data(), isft_length);

	// sdta
	file.read32(); // LIST
	file.read32(); // LIST chunk length
	file.read32(); // sdta
	file.read32(); // smpl
	uint32_t smpl_length = file.read32();

	std::vector<uint8_t> raw_samples = file.readn(smpl_length, true);
	int range = pow(2, 16);
	float max = range / 2.0 - 1.0;

	for (int i = 0; i < raw_samples.size(); i += 2) {
		uint8_t byte1 = raw_samples[i];
		uint8_t byte2 = raw_samples[i + 1];

		int16_t sample = byte2 << 8 | byte1;

		samples.push_back(sample / max);
	}

	std::ofstream audio("audio.bin", std::ios::binary);

	for (uint8_t smpl : raw_samples) {
		audio << smpl;
	}

	audio.close();

	// pdta
	file.read32();                        // LIST
	file.read32();                        // LIST chunk length
	file.read32();                        // pdta
	file.read32();                        // phdr
	uint32_t phdr_length = file.read32(); // phdr chunk length

	for (int i = 0; i < phdr_length; i += 38) {
		sf_preset_header preset_header;
		std::vector<uint8_t> preset_name = file.readn(20, true);
		memcpy(preset_header.preset_name, preset_name.data(), 20);
		preset_header.preset = file.read16();
		preset_header.bank = file.read16();
		preset_header.perset_bag_ndx = file.read16();
		preset_header.library = file.read32();
		preset_header.genre = file.read32();
		preset_header.morphology = file.read32();
		preset_headers.push_back(preset_header);
	}

	file.read32();                        // pbag
	uint32_t pbag_length = file.read32(); // pbag chunk length
	for (int i = 0; i < pbag_length; i += 4) {
		sf_preset_bag preset_index;
		preset_index.GenNdx = file.read16();
		preset_index.GenNdx = file.read16();
		preset_index_list.push_back(preset_index);
	}
	// file.add_to_cursor(pbag_length);

	file.read32();                        // pmod
	uint32_t pmod_length = file.read32(); // pmod chunk length
	for (int i = 0; i < pmod_length; i += 10) {
		sf_preset_mod_list preset_modulator;

		preset_modulator.src_oper = (sf_modulator)file.read16();
		preset_modulator.dest_oper = (sf_generator)file.read16();
		preset_modulator.mod_amount = file.read16();
		preset_modulator.amt_src_oper = (sf_modulator)file.read16();
		preset_modulator.trans_oper = (sf_transform)file.read16();

		preset_modulator_list.push_back(preset_modulator);
	}
	// file.add_to_cursor(pmod_length);

	file.read32();                        // pgen
	uint32_t pgen_length = file.read32(); // pgen chunk length
	for (int i = 0; i < pgen_length; i += 4) {
		sf_preset_gen_list preset_generator;

		preset_generator.oper = (sf_generator)file.read16();
		preset_generator.amount.w_amount = file.read16();

		preset_generate_list.push_back(preset_generator);
	}
	// file.add_to_cursor(pgen_length);

	file.read32();                        // inst
	uint32_t inst_length = file.read32(); // inst chunk length
	for (int i = 0; i < inst_length; i += 22) {
		sf_inst instument;

		memcpy(instument.inst_name, file.readn(20, true).data(), 20);
		instument.inst_bag_ndx = file.read16();

		instrument_names_indices.push_back(instument);
	}
	// file.add_to_cursor(inst_length);

	file.read32();                        // ibag
	uint32_t ibag_length = file.read32(); // ibag chunk length
	for (int i = 0; i < ibag_length; i += 4) {
		sf_inst_bag instument_bag;

		instument_bag.inst_gen_ndx = file.read16();
		instument_bag.inst_mod_ndx = file.read16();

		instument_index_list.push_back(instument_bag);
	}
	// file.add_to_cursor(ibag_length);

	file.read32();                        // imod
	uint32_t imod_length = file.read32(); // imod chunk length
	for (int i = 0; i < imod_length; i += 10) {
		sf_inst_mod_list instrument_modulator;

		instrument_modulator.src_oper = (sf_modulator)file.read16();
		instrument_modulator.dest_oper = (sf_generator)file.read16();
		instrument_modulator.mod_amount = file.read16();
		instrument_modulator.amt_src_oper = (sf_modulator)file.read16();
		instrument_modulator.trans_oper = (sf_transform)file.read16();

		instrument_modulator_list.push_back(instrument_modulator);
	}
	// file.add_to_cursor(imod_length);

	file.read32();                        // igen
	uint32_t igen_length = file.read32(); // igen chunk length
	for (int i = 0; i < igen_length; i += 4) {
		sf_inst_gen_list instrument_generator;

		instrument_generator.oper = (sf_generator)file.read16();
		instrument_generator.amount.w_amount = file.read16();

		instrument_generator_list.push_back(instrument_generator);
	}
	// file.add_to_cursor(igen_length);

	file.read32();                        // shdr
	uint32_t shdr_length = file.read32(); // shdr chunk length
	for (int i = 0; i < shdr_length; i += 46) {
		sf_sample sample_header;
		std::vector<uint8_t> sample_name = file.readn(20, true);
		memcpy(sample_header.sample_name, sample_name.data(), 20);
		sample_header.start = file.read32();
		sample_header.end = file.read32();
		sample_header.start_loop = file.read32();
		sample_header.end_loop = file.read32();
		sample_header.sample_rate = file.read32();
		sample_header.original_key = file.read8();
		sample_header.correction = file.read8();
		sample_header.sample_link = file.read16();
		std::cout << "test " << shdr_length / 46 << "/" << (i / 46) + 1 << "\n";
		sample_header.sample_type = (sf_sample_link)file.read16();

		sample_headers.push_back(sample_header);
	}

	// std::cout << std::hex << file.read32() << "\n";
}
