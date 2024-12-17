#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <ios>
#include <iostream>
#include <optional>
#include <sys/types.h>
#include <variant>
#include <vector>

std::ofstream log_file("./log.txt");

struct MIDI_DIVISION {
    bool format;
    union {
        uint16_t ticks_per_quarter_note;
        struct {
            uint8_t frames_per_second;
            uint8_t ticks_per_frame;
        };
    };
};

struct MIDI_HEADER {
    char type[5] = {0};
    uint32_t length;
    uint16_t format;
    uint16_t ntracks;
    MIDI_DIVISION division;
};

struct MIDI_SYSEX_EVENT {
    std::vector<uint8_t> data;
};

struct MIDI_META_EVENT {
    uint8_t type;
    std::vector<uint8_t> data;
};

struct MIDI_MIDI_EVENT {
    uint8_t type;
    uint8_t channel;
    std::vector<uint8_t> data;
};

struct MIDI_EVENT {
    uint32_t delta_time;
    uint8_t type;
    std::variant<MIDI_SYSEX_EVENT, MIDI_META_EVENT, MIDI_MIDI_EVENT> event;
};

struct MIDI_TRACK {
    char type[5] = {0};
    uint32_t length;
    std::vector<MIDI_EVENT> events;
};

struct MIDI_FILE {
    MIDI_HEADER header;
    std::vector<MIDI_TRACK> tracks;
};

uint8_t read_8(std::ifstream *file) {
    char byte;

    file->get(byte);

    return byte;
}

uint16_t read_16(std::ifstream *file) {
    char bytes[2];

    file->read(bytes, 2);

    uint8_t byte1 = bytes[0];
    uint8_t byte2 = bytes[1];

    return byte1 << 8 | byte2;
}

uint32_t read_32(std::ifstream *file) {
    char bytes[4];

    file->read(bytes, 4);

    uint8_t byte1 = bytes[0];
    uint8_t byte2 = bytes[1];
    uint8_t byte3 = bytes[2];
    uint8_t byte4 = bytes[3];

    return byte1 << 24 | byte2 << 16 | byte3 << 8 | byte4;
}

uint32_t read_number(std::ifstream *file) {
    uint32_t number = 0;
    char byte = 0;
    bool end = false;
    file->get(byte);

    while (!end) {

        if ((byte & 0b10000000) > 0) {
            number = (number << 7) + (byte & 0b1111111);
            file->get(byte);
        } else {
            number = (number << 7) + byte;
            end = true;
        }
    }

    return number;
}

MIDI_MIDI_EVENT read_midi_midi_event(std::ifstream *file, uint8_t typechannel) {
    MIDI_MIDI_EVENT event;

    event.type = typechannel & 0xF0;
    event.channel = typechannel & 0x0F;

    switch (event.type >> 4) {
    case 0x8:
    case 0x9: {
        uint32_t key_number = read_8(file);
        uint32_t velocity = read_8(file);
        event.data.push_back(key_number);
        event.data.push_back(velocity);
    } break;
    case 0xB: {
        uint8_t controller_number = read_8(file);
        uint8_t value = read_8(file);
        event.data.push_back(controller_number);
        event.data.push_back(value);
    } break;
    case 0xC: {
        uint8_t program_number = read_8(file);
        event.data.push_back(program_number);
    } break;
    }

    return event;
};

std::vector<MIDI_EVENT> read_midi_track_events(std::ifstream *file,
                                               uint32_t end_location) {
    std::vector<MIDI_EVENT> events;

    while (file->tellg() < end_location) {
        uint32_t delta_time = read_number(file);
        uint8_t event_type = read_8(file);
        MIDI_EVENT event;
        switch (event_type) {
        case 0xF0: {
            MIDI_SYSEX_EVENT sysex_event;
            uint32_t length = read_number(file);
            for (uint32_t i = 0; i < length; i++) {
                sysex_event.data.push_back(read_8(file));
            }
            event.event = sysex_event;
        } break;
        case 0xFF: {
            MIDI_META_EVENT meta_event;
            uint8_t meta_type = read_8(file);
            uint32_t event_length = read_number(file);
            uint32_t end = (uint32_t)file->tellg() + event_length;
            meta_event.type = meta_type;
            for (uint32_t i = 0; i < event_length; i++) {
                uint8_t byte = read_8(file);
                meta_event.data.push_back(byte);
            }
            event.event = meta_event;

            file->seekg(end);
        } break;
        case 0x80 ... 0xEF: {
            event.event = read_midi_midi_event(file, event_type);
        } break;
        default: {
            log_file << "delta: " << delta_time << ", type:" << std::hex
                     << (uint)event_type << std::dec << "("
                     << std::bitset<8>(event_type) << ")" << "\n";
            return events; // exits loop
        }
        }
        events.push_back(event);
    }

    return events;
}

std::optional<MIDI_FILE> read_midi_file(const char *filename) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (file.fail()) {
        std::cerr << "couldn't open file " << filename << "\n";
        return std::nullopt;
    }
    MIDI_FILE midi_file = {};
    file.read(midi_file.header.type, 4);
    midi_file.header.length = read_32(&file);
    midi_file.header.format = read_16(&file);
    midi_file.header.ntracks = read_16(&file);
    uint16_t division_number = read_16(&file);
    std::bitset<16> division = division_number;
    bool division_format = division[15];
    midi_file.header.division.format = division_format;
    if (division_format) {
        midi_file.header.division.ticks_per_frame = division_number & 0x00FF;
        midi_file.header.division.frames_per_second =
            256 - ((division_number >> 8) &
                   0x00FF); // gets absolute of the number. -30 -> 30
    } else {
        midi_file.header.division.ticks_per_quarter_note = division_number;
    }
    midi_file.tracks = std::vector<MIDI_TRACK>();
    midi_file.tracks.reserve(midi_file.header.ntracks);

    for (uint16_t i = 0; i < midi_file.header.ntracks; i++) {
        MIDI_TRACK track = {};
        file.read(track.type, 4);
        track.length = read_32(&file);
        uint32_t end_location = (uint32_t)file.tellg() + track.length;
        track.events = read_midi_track_events(&file, end_location);
        file.seekg(end_location);
        midi_file.tracks.push_back(track);
    }

    return midi_file;
}

void print_event(MIDI_MIDI_EVENT event) {
    std::cout << "MIDI:\n"
              << " type: 0x" << std::hex << std::uppercase
              << (uint)(event.type >> 4) << std::dec << std::nouppercase << "("
              << std::bitset<4>(event.type >> 4) << ")"
              << "\n channel: " << (uint)event.channel << "\n data: \n";
    switch (event.type >> 4) {
    case 0x8:
    case 0x9: {
        uint32_t key_number = event.data[0];
        uint32_t velocity = event.data[1];
        if (event.type == 0x80 || (event.type == 0x90 && velocity == 0)) {
            std::cout << "  note: OFF\n";
        } else {
            std::cout << "  note: ON\n";
        }
        std::cout << "  key number: " << key_number
                  << "\n  velocity: " << velocity << "\n";
    } break;
    case 0xB: {
        uint8_t controller_number = event.data[0];
        uint8_t value = event.data[1];
        std::cout << "  controller number: " << (uint)controller_number
                  << "\n  value: " << (uint)value << "\n";

    } break;
    case 0xC: {
        uint8_t program_number = event.data[0];
        std::cout << "  program change: " << (uint)program_number << "\n";
    } break;
    }
}

void print_event(MIDI_META_EVENT event) {
    std::cout << "META:\n"
              << " meta type: 0x" << std::hex << std::uppercase
              << (uint)event.type << std::dec << std::nouppercase
              << "\n length: " << event.data.size() << "\n data:\n"
              << std::hex << std::uppercase << "  ";

    for (uint8_t byte : event.data) {
        if (event.type >= 0x01 && event.type <= 0x07) {
            std::cout << byte;
        } else {
            std::cout << (uint)byte << " ";
        }
    }
    std::cout << std::dec << std::nouppercase << "\n";
}

void print_event(MIDI_SYSEX_EVENT event) {
    std::cout << "SYSEX:\n data:\n  " << std::hex << std::uppercase;

    for (uint8_t byte : event.data) {
        std::cout << (uint)byte << " ";
    }

    std::cout << std::dec << std::nouppercase << "\n";
}

const char *filename = "./sample/test2.mid";

int main(int argc, char **argv) {
    std::cout << "current file:" << filename << "\n";

    std::optional<MIDI_FILE> file_optional = read_midi_file(filename);

    if (!file_optional.has_value()) {
        return 1;
    }

    MIDI_FILE file = file_optional.value();
    std::cout << "header: " << file.header.type << "\n";
    std::cout << "header length: " << file.header.length << "\n";
    std::cout << "format: " << file.header.format << "\n";
    std::cout << "tracks: " << file.header.ntracks << "\n";
    std::cout << "dividion type: " << file.header.division.format << "\n";
    if (file.header.division.format) {
        std::cout << "dividion ticks per frame: "
                  << file.header.division.ticks_per_frame << "\n";
        std::cout << "dividion frames per second: "
                  << file.header.division.frames_per_second << "\n";
    } else {
        std::cout << "dividion ticks per quarter note: "
                  << file.header.division.ticks_per_quarter_note << "\n";
    }
    std::cout << "------------------------------------------------------\n";
    for (uint16_t i = 0; i < file.header.ntracks; i++) {
        std::cout << "track " << i << "\n";
        MIDI_TRACK track = file.tracks[i];

        for (MIDI_EVENT midi_event : track.events) {
            std::cout << "\n";

            std::cout << "delta: " << (uint)midi_event.delta_time << "\n";

            if (std::holds_alternative<MIDI_MIDI_EVENT>(midi_event.event)) {
                MIDI_MIDI_EVENT event =
                    std::get<MIDI_MIDI_EVENT>(midi_event.event);
                print_event(event);
            } else if (std::holds_alternative<MIDI_META_EVENT>(
                           midi_event.event)) {
                MIDI_META_EVENT event =
                    std::get<MIDI_META_EVENT>(midi_event.event);
                print_event(event);
            } else if (std::holds_alternative<MIDI_SYSEX_EVENT>(
                           midi_event.event)) {
                MIDI_SYSEX_EVENT event =
                    std::get<MIDI_SYSEX_EVENT>(midi_event.event);
                print_event(event);
            }
        }

        std::cout << "------------------------------------------------------\n";
    }

    return 0;
}
