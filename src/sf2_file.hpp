#ifndef SF2_FILE_HPP_
#define SF2_FILE_HPP_
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cwchar>
#include <ostream>
#include <random>
#include <string>
#include <vector>

struct sf_ranges_type {
	uint8_t low;
	uint8_t high;
};

union sf_gen_amount_type {
	sf_ranges_type ranges;
	int16_t sh_amount;
	uint16_t w_amount;
};

enum sf_sample_link {
	mono_sample = 1,
	right_sample = 2,
	left_sample = 4,
	linked_sample = 8,
	rom_mono_sample = 0x8001,
	rom_right_sample = 0x8002,
	rom_left_sample = 0x8004,
	rom_linked_sample = 0x8008,
};

struct sf_preset_header { // sfPresetHeader
	char preset_name[20];
	uint16_t preset;
	uint16_t bank;
	uint16_t perset_bag_ndx;
	uint32_t library;
	uint32_t genre;
	uint32_t morphology;
};

struct sf_preset_bag { // sfPresetBag
	uint16_t GenNdx;
	uint16_t ModNdx;
};

typedef uint16_t sf_modulator; // enum

enum class sf_generator : uint16_t {
	startAddrOffset = 0,
	endAddrOffset = 1,
	startloopAddrOffset = 2,
	endloopAddrOffset = 3,
	startAddrCoarseOffset = 4,
	modLfoToPitch = 5,
	vibLfoToPitch = 6,
	modEnvToPitch = 7,
	initialFilterFc = 8,
	initialFilterQ = 9,
	modLfoToFilterFc = 10,
	modEnvToFilterFc = 11,
	endAddrCoarseOffset = 12,
	modLfoToVolume = 13,
	chorusEffectsSend = 15,
	reverbEffectsSend = 16,
	pan = 17,
	delayModLFO = 21,
	freqModLFO = 22,
	delayVibLFO = 23,
	freqVibLFO = 24,
	delayModEnv = 25,
	attackModEnv = 26,
	holdModEnv = 27,
	decayModEnv = 28,
	sustainModEnv = 29,
	releaseModEnv = 30,
	keynumToModEnvHold = 31,
	keynumToModEnvDecay = 32,
	delayVolEnv = 33,
	attackVolEnv = 34,
	holdVolEnv = 35,
	decayVolEnv = 36,
	sustainVolEnv = 37,
	releaseVolEnv = 38,
	keynumToVolEnvHold = 39,
	keynumToVolEnvDecay = 40,
	instrument = 41,
	keyRange = 43,
	velRange = 44,
	startloopAddrCoarseOffset = 45,
	keynum = 46,
	velocity = 47,
	initialAttenuation = 48,
	endloopAddrCoarseOffset = 50,
	coarseTune = 51,
	fineTune = 52,
	sampleID = 53,
	sampleModes = 54,
	scaleTuning = 56,
	exclusiveClass = 57,
	overridingRootKey = 58,
	endOper = 60
};

std::ostream &operator<<(std::ostream &stream, const sf_generator &gen);

enum class sf_transform : uint16_t {
	linear = 0,
	absolute = 2,
};

struct sf_preset_mod { // sfModList
	sf_modulator src_oper;
	sf_generator dest_oper;
	int16_t mod_amount;
	sf_modulator amt_src_oper;
	sf_transform trans_oper;
};

struct sf_preset_gen {
	sf_generator oper;
	sf_gen_amount_type amount;
};

struct sf_inst {
	char inst_name[20];
	uint16_t inst_bag_ndx;
};

struct sf_inst_bag {
	uint16_t inst_gen_ndx;
	uint16_t inst_mod_ndx;
};

typedef sf_preset_mod sf_inst_mod;
typedef sf_preset_gen sf_inst_gen;

struct sf_sample {
	char sample_name[20];
	uint32_t start;
	uint32_t end;
	uint32_t start_loop;
	uint32_t end_loop;
	uint32_t sample_rate;
	uint8_t original_key;
	char correction;
	uint16_t sample_link;
	sf_sample_link sample_type;
};

class SF2File {
  public:
	SF2File(std::string file_path);

  private:
  public:
	// INFO
	uint32_t major_version;
	uint32_t minor_version;
	std::string sound_engine;
	std::string name;
	std::string date_of_creation;
	std::string tool;
	// sdta
	std::vector<float> samples;
	// pdta
	std::vector<sf_preset_header> preset_headers;
	std::vector<sf_preset_bag> preset_index_list;
	std::vector<sf_preset_mod> preset_modulator_list;
	std::vector<sf_preset_gen> preset_generator_list;
	std::vector<sf_inst> instrument_names_indices;
	std::vector<sf_inst_bag> instument_index_list;
	std::vector<sf_inst_mod> instrument_modulator_list;
	std::vector<sf_inst_gen> instrument_generator_list;
	std::vector<sf_sample> sample_headers;

  private:
	std::string file_path;
};

#endif // SF2_FILE_HPP
