#include "aare/Dtype.hpp"
#include "aare/logger.hpp"
#include "helpers/FileInterface.hpp"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace aare {

class CustomMythenFile : public DetectorFileInterface {

  public:
    CustomMythenFile(const ssize_t rows, const ssize_t cols,
                     const std::string &mode = "r");

    ~CustomMythenFile();

    Frame read_frame() override;

    void read_into(std::byte *image_buf) override;

    size_t bytes_per_frame() override;

    size_t pixels_per_frame() override;

  private:
    // uint8_t m_num_counts{}; TODO extend
    static const Dtype m_dtype;
    static const DetectorType m_det_type = DetectorType::Mythen3;
};

class CustomBadChannelsFile {

  public:
    CustomBadChannelsFile(const std::string &filename) : m_filename(filename) {
        try {
            m_file = std::ifstream(filename, std::ios_base::in);
            if (!m_file.good()) {
                throw std::logic_error("file does not exist");
            }
        } catch (std::exception &e) {
            LOG(TLogLevel::logERROR) << "file does not exist\n";
        }
    }

    ~CustomBadChannelsFile() { m_file.close(); }

    void read_into_array(NDView<bool, 1> array) {
        std::string line;
        try {
            while (std::getline(m_file, line)) {
                std::size_t pos = line.find("-");

                if (pos == std::string::npos) {
                    array(std::stoi(line)) = true;
                } else {
                    size_t line_size = line.size();
                    for (int i = std::stoi(line.substr(0, pos));
                         i <= std::stoi(line.substr(pos + 1, line_size - pos));
                         ++i)
                        array(i) = true;
                }
            }
        } catch (const std::exception &e) {
            LOG(TLogLevel::logERROR) << e.what();
        }
    }

  private:
    std::string m_filename{};
    std::ifstream m_file{};
};

} // namespace aare