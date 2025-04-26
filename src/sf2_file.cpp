#include "sf2_file.hpp"
#include "file.hpp"
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <ios>
#include <iostream>
#include <string>
#include <vector>

void readListChunk(SF2File *sf2_file, File &file);
void readINFOChunk(SF2File *sf2_file, File &file);
void readSDTAChunk(SF2File *sf2_file, File &file);
void readPDTAChunk(SF2File *sf2_file, File &file);

SF2File::SF2File(std::string file_name) {
	file_path = file_name;
	File file(file_name);

	file.read32(); // RIFF
	file.read32(); // riff chunk length
	file.read32(); // SFBK

	while (file.get_cursor() < file.get_file_length()) {
		readListChunk(this, file);
	}
}

void readListChunk(SF2File *sf2_file, File &file) {
	file.read32();                   // LIST
	uint32_t length = file.read32(); // LIST length
	uint32_t end_of_chunk = file.get_cursor() + length;
	uint32_t type = file.read32();
	uint8_t *typechars = (uint8_t *)&type;

	switch (type) {
	case 0x4F464E49: // INFO
		readINFOChunk(sf2_file, file);
		break;
	case 0x61746473: // sdta
		readSDTAChunk(sf2_file, file);
		break;
	case 0x61746470: // pdta
		readPDTAChunk(sf2_file, file);
		break;
	}

	file.set_cursor(end_of_chunk);
}

void readINFOChunk(SF2File *sf2_file, File &file) {
	file.read32(); // ifil
	file.read32(); // ifil length
	sf2_file->major_version = file.read16();
	sf2_file->minor_version = file.read16();
	file.read32();                        // isng
	uint32_t isng_length = file.read32(); // isng length
	sf2_file->sound_engine.assign((char *)file.readn(isng_length, true).data(),
	                              isng_length);
	file.read32(); // INAM
	uint32_t inam_length = file.read32();
	sf2_file->name.assign((char *)file.readn(inam_length, true).data(),
	                      inam_length);
	file.read32(); // ICRD
	uint32_t icrd_length = file.read32();
	sf2_file->date_of_creation.assign(
	    (char *)file.readn(icrd_length, true).data(), icrd_length);
	file.read32(); // ISFT
	uint32_t isft_length = file.read32();
	sf2_file->tool.assign((char *)file.readn(isft_length, true).data(),
	                      isft_length);
}

void readSDTAChunk(SF2File *sf2_file, File &file) {
	file.read32(); // smpl
	uint32_t smpl_length = file.read32();

	std::vector<uint8_t> raw_samples = file.readn(smpl_length, true);
	int range = pow(2, 16);
	float max = range / 2.0 - 1.0;

	for (int i = 0; i < raw_samples.size(); i += 2) {
		uint8_t byte1 = raw_samples[i];
		uint8_t byte2 = raw_samples[i + 1];

		int16_t sample = byte2 << 8 | byte1;

		sf2_file->samples.push_back(sample / max);
	}

	std::ofstream audio("audio.bin", std::ios::binary);

	for (uint8_t smpl : raw_samples) {
		audio << smpl;
	}

	audio.close();
}

void readPDTAChunk(SF2File *sf2_file, File &file) {
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
		sf2_file->preset_headers.push_back(preset_header);
	}

	file.read32();                        // pbag
	uint32_t pbag_length = file.read32(); // pbag chunk length
	for (int i = 0; i < pbag_length; i += 4) {
		sf_preset_bag preset_index;
		preset_index.GenNdx = file.read16();
		preset_index.ModNdx = file.read16();
		sf2_file->preset_index_list.push_back(preset_index);
	}
	// file.add_to_cursor(pbag_length);

	file.read32();                        // pmod
	uint32_t pmod_length = file.read32(); // pmod chunk length
	for (int i = 0; i < pmod_length; i += 10) {
		sf_preset_mod preset_modulator;

		preset_modulator.src_oper = (sf_modulator)file.read16();
		preset_modulator.dest_oper = (sf_generator)file.read16();
		preset_modulator.mod_amount = file.read16();
		preset_modulator.amt_src_oper = (sf_modulator)file.read16();
		preset_modulator.trans_oper = (sf_transform)file.read16();

		sf2_file->preset_modulator_list.push_back(preset_modulator);
	}
	// file.add_to_cursor(pmod_length);

	file.read32();                        // pgen
	uint32_t pgen_length = file.read32(); // pgen chunk length
	for (int i = 0; i < pgen_length; i += 4) {
		sf_preset_gen preset_generator;

		preset_generator.oper = (sf_generator)file.read16();
		preset_generator.amount.w_amount = file.read16();

		sf2_file->preset_generator_list.push_back(preset_generator);
	}
	// file.add_to_cursor(pgen_length);

	file.read32();                        // inst
	uint32_t inst_length = file.read32(); // inst chunk length
	for (int i = 0; i < inst_length; i += 22) {
		sf_inst instument;

		memcpy(instument.inst_name, file.readn(20, true).data(), 20);
		instument.inst_bag_ndx = file.read16();

		sf2_file->instrument_names_indices.push_back(instument);
	}
	// file.add_to_cursor(inst_length);

	file.read32();                        // ibag
	uint32_t ibag_length = file.read32(); // ibag chunk length
	for (int i = 0; i < ibag_length; i += 4) {
		sf_inst_bag instument_bag;

		instument_bag.inst_gen_ndx = file.read16();
		instument_bag.inst_mod_ndx = file.read16();

		sf2_file->instument_index_list.push_back(instument_bag);
	}
	// file.add_to_cursor(ibag_length);

	file.read32();                        // imod
	uint32_t imod_length = file.read32(); // imod chunk length
	for (int i = 0; i < imod_length; i += 10) {
		sf_inst_mod instrument_modulator;

		instrument_modulator.src_oper = (sf_modulator)file.read16();
		instrument_modulator.dest_oper = (sf_generator)file.read16();
		instrument_modulator.mod_amount = file.read16();
		instrument_modulator.amt_src_oper = (sf_modulator)file.read16();
		instrument_modulator.trans_oper = (sf_transform)file.read16();

		sf2_file->instrument_modulator_list.push_back(instrument_modulator);
	}
	// file.add_to_cursor(imod_length);

	file.read32();                        // igen
	uint32_t igen_length = file.read32(); // igen chunk length
	for (int i = 0; i < igen_length; i += 4) {
		sf_inst_gen instrument_generator;

		instrument_generator.oper = (sf_generator)file.read16();
		instrument_generator.amount.w_amount = file.read16();

		sf2_file->instrument_generator_list.push_back(instrument_generator);
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
		sample_header.sample_type = (sf_sample_link)file.read16();

		sf2_file->sample_headers.push_back(sample_header);
	}
}

std::ostream &operator<<(std::ostream &stream, const sf_generator &gen) {
	switch (gen) {
	case sf_generator::startAddrOffset:
		stream << "startAddrOffset";
		break;
	case sf_generator::endAddrOffset:
		stream << "endAddrOffset";
		break;
	case sf_generator::startloopAddrOffset:
		stream << "startloopAddrOffset";
		break;
	case sf_generator::endloopAddrOffset:
		stream << "endloopAddrOffset";
		break;
	case sf_generator::startAddrCoarseOffset:
		stream << "startAddrCoarseOffset";
		break;
	case sf_generator::modLfoToPitch:
		stream << "modLfoToPitch";
		break;
	case sf_generator::vibLfoToPitch:
		stream << "vibLfoToPitch";
		break;
	case sf_generator::modEnvToPitch:
		stream << "modEnvToPitch";
		break;
	case sf_generator::initialFilterFc:
		stream << "initialFilterFc";
		break;
	case sf_generator::initialFilterQ:
		stream << "initialFilterQ";
		break;
	case sf_generator::modLfoToFilterFc:
		stream << "modLfoToFilterFc";
		break;
	case sf_generator::modEnvToFilterFc:
		stream << "modEnvToFilterFc";
		break;
	case sf_generator::endAddrCoarseOffset:
		stream << "endAddrCoarseOffset";
		break;
	case sf_generator::modLfoToVolume:
		stream << "modLfoToVolume";
		break;
	case sf_generator::chorusEffectsSend:
		stream << "chorusEffectsSend";
		break;
	case sf_generator::reverbEffectsSend:
		stream << "reverbEffectsSend";
		break;
	case sf_generator::pan:
		stream << "pan";
		break;
	case sf_generator::delayModLFO:
		stream << "delayModLFO";
		break;
	case sf_generator::freqModLFO:
		stream << "freqModLFO";
		break;
	case sf_generator::delayVibLFO:
		stream << "delayVibLFO";
		break;
	case sf_generator::freqVibLFO:
		stream << "freqVibLFO";
		break;
	case sf_generator::delayModEnv:
		stream << "delayModEnv";
		break;
	case sf_generator::attackModEnv:
		stream << "attackModEnv";
		break;
	case sf_generator::holdModEnv:
		stream << "holdModEnv";
		break;
	case sf_generator::decayModEnv:
		stream << "decayModEnv";
		break;
	case sf_generator::sustainModEnv:
		stream << "sustainModEnv";
		break;
	case sf_generator::releaseModEnv:
		stream << "releaseModEnv";
		break;
	case sf_generator::keynumToModEnvHold:
		stream << "keynumToModEnvHold";
		break;
	case sf_generator::keynumToModEnvDecay:
		stream << "keynumToModEnvDecay";
		break;
	case sf_generator::delayVolEnv:
		stream << "delayVolEnv";
		break;
	case sf_generator::attackVolEnv:
		stream << "attackVolEnv";
		break;
	case sf_generator::holdVolEnv:
		stream << "holdVolEnv";
		break;
	case sf_generator::decayVolEnv:
		stream << "decayVolEnv";
		break;
	case sf_generator::sustainVolEnv:
		stream << "sustainVolEnv";
		break;
	case sf_generator::releaseVolEnv:
		stream << "releaseVolEnv";
		break;
	case sf_generator::keynumToVolEnvHold:
		stream << "keynumToVolEnvHold";
		break;
	case sf_generator::keynumToVolEnvDecay:
		stream << "keynumToVolEnvDecay";
		break;
	case sf_generator::instrument:
		stream << "instrument";
		break;
	case sf_generator::keyRange:
		stream << "keyRange";
		break;
	case sf_generator::velRange:
		stream << "velRange";
		break;
	case sf_generator::startloopAddrCoarseOffset:
		stream << "startloopAddrCoarseOffset";
		break;
	case sf_generator::keynum:
		stream << "keynum";
		break;
	case sf_generator::velocity:
		stream << "velocity";
		break;
	case sf_generator::initialAttenuation:
		stream << "initialAttenuation";
		break;
	case sf_generator::endloopAddrCoarseOffset:
		stream << "endloopAddrCoarseOffset";
		break;
	case sf_generator::coarseTune:
		stream << "coarseTune";
		break;
	case sf_generator::fineTune:
		stream << "fineTune";
		break;
	case sf_generator::sampleID:
		stream << "sampleID";
		break;
	case sf_generator::sampleModes:
		stream << "sampleModes";
		break;
	case sf_generator::scaleTuning:
		stream << "scaleTuning";
		break;
	case sf_generator::exclusiveClass:
		stream << "exclusiveClass";
		break;
	case sf_generator::overridingRootKey:
		stream << "overridingRootKey";
		break;
	case sf_generator::endOper:
		stream << "endOper";
		break;
	}
	return stream;
}
