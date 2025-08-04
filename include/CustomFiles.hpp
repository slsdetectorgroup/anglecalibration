#pragma once

#include "aare/Dtype.hpp"
#include "helpers/FileInterface.hpp"
#include "logger.hpp"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace angcal {

class CustomMythenFile : public SimpleFileInterface {
  public:
    CustomMythenFile() = default;

    ~CustomMythenFile() { m_file.close(); }

    void read_into(std::byte *image_buf,
                   const ssize_t data_types_bytes = 4) override {
        uint32_t strip_index, photon_count;
        try {
            while (m_file >> strip_index >> photon_count) {
                std::memcpy(image_buf, &photon_count, sizeof(photon_count));

                image_buf += data_types_bytes;
            }
        } catch (const std::exception &e) {
            LOG(TLogLevel::logERROR) << e.what();
        }
    }
};

class CustomFlatFieldFile : public SimpleFileInterface {
  public:
    CustomFlatFieldFile() = default;

    ~CustomFlatFieldFile() { m_file.close(); }

    void read_into(std::byte *image_buf,
                   const ssize_t data_types_bytes = 8) override {
        uint32_t strip_index;
        double flatfield_value, flatfield_error, val_1, val_2;
        std::string line;

        try {
            std::getline(m_file, line); // skip header
            while (m_file >> strip_index >> flatfield_value >>
                   flatfield_error >> val_1 >> val_2) {
                std::memcpy(image_buf, &flatfield_value,
                            sizeof(flatfield_value));

                image_buf += data_types_bytes;
            }
        } catch (const std::exception &e) {
            LOG(TLogLevel::logERROR) << e.what();
        }
    }
};

class CustomBadChannelsFile : public SimpleFileInterface {

  public:
    CustomBadChannelsFile() = default;

    ~CustomBadChannelsFile() { m_file.close(); }

    void read_into(std::byte *buffer, const ssize_t data_type_bytes = 1) {
        std::string line;

        size_t index = 0;
        try {
            while (std::getline(m_file, line)) {

                std::size_t pos = line.find("-");

                if (pos == std::string::npos) {
                    index = std::stoi(line) * data_type_bytes;
                    buffer[index] = static_cast<std::byte>(true);
                } else {
                    size_t line_size = line.size();
                    for (int i = std::stoi(line.substr(0, pos));
                         i <= std::stoi(line.substr(pos + 1, line_size - pos));
                         ++i) {
                        index = i * data_type_bytes;
                        buffer[index] = static_cast<std::byte>(true);
                    }
                }
            }
        } catch (const std::exception &e) {
            LOG(TLogLevel::logERROR) << e.what();
        }
    }
};

class InitialAngCalParametersFile : public SimpleFileInterface {
  public:
    InitialAngCalParametersFile() = default;
    ~InitialAngCalParametersFile() { m_file.close(); }

    void read_into(std::byte *buffer, const ssize_t data_type_bytes = 8) {
        std::string line;
        uint32_t module_number{};

        try {
            size_t index = 0;

            while (m_file >> line) {
                if (line == "module") {
                    m_file >> line;
                    module_number = std::stoi(line);
                    index = module_number * data_type_bytes * 3;
                }
                if (line == "center") {
                    m_file >> line;

                    *reinterpret_cast<double *>(&buffer[index]) =
                        std::stod(line);
                }
                if (line == "conversion") {
                    m_file >> line;
                    *reinterpret_cast<double *>(
                        &buffer[index + data_type_bytes]) = std::stod(line);
                }
                if (line == "offset") {
                    m_file >> line;
                    *reinterpret_cast<double *>(
                        &buffer[index + 2 * data_type_bytes]) = std::stod(line);
                }
            }
        } catch (const std::exception &e) {
            LOG(TLogLevel::logERROR) << "Error: " << e.what() << std::endl;
        }
    }
};

} // namespace angcal