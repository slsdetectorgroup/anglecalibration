#include "aare/Dtype.hpp"
#include "aare/logger.hpp"
#include "helpers/FileInterface.hpp"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace angcal {

class CustomMythenFile : public DetectorFileInterface {

  public:
    CustomMythenFile(const ssize_t rows, const ssize_t cols,
                     const std::string &mode = "r");

    ~CustomMythenFile();

    aare::Frame read_frame() override;

    void read_into(std::byte *image_buf) override;

    size_t bytes_per_frame() override;

    size_t pixels_per_frame() override;

  private:
    // uint8_t m_num_counts{}; TODO extend
    static const aare::Dtype m_dtype;
    static const aare::DetectorType m_det_type = aare::DetectorType::Mythen3;
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
            LOG(aare::TLogLevel::logERROR) << e.what();
        }
    }
};

} // namespace angcal