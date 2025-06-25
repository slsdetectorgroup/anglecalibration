#include "aare/Dtype.hpp"
#include "aare/FileInterface.hpp"
#include "aare/logger.hpp"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace aare {

// TODO a lot to overload
class CustomMythenFile : public FileInterface {

  public:
    CustomMythenFile(const std::string &filename, const ssize_t rows,
                     const ssize_t cols = 1, const std::string &mode = "r");

    ~CustomMythenFile();

    Frame read_frame() override;

    Frame read_frame(size_t frame_number) override;

    std::vector<Frame> read_n(size_t n_frames) override;

    void read_into(std::byte *image_buf) override;

    void read_into(std::byte *image_buf, size_t n_frames) override;

    size_t frame_number(size_t frame_index) override;

    size_t bytes_per_frame() override;

    size_t pixels_per_frame() override;

    void seek(size_t frame_number) override;

    size_t tell() override;

    size_t total_frames() const override;

    size_t rows() const override;
    size_t cols() const override;

    size_t bitdepth() const override;

    DetectorType detector_type() const override;

  private:
    std::string m_filename{};
    std::ifstream m_file{};
    ssize_t m_num_strips{};
    // uint8_t m_num_counts{}; TODO extend
    ssize_t m_rows{};
    ssize_t m_cols{};
    static const Dtype m_dtype;
    static const DetectorType det_type = DetectorType::Mythen3;
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