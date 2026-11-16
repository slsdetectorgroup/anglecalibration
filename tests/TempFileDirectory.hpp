#pragma once
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

/**
 * @brief creates a unique temporary directory which is removed (including
 * all files written into it) on destruction
 */
class TemporaryDirectory {
  public:
    explicit TemporaryDirectory(const std::string &prefix = "angcal-utils") {
        const auto unique =
            std::chrono::steady_clock::now().time_since_epoch().count();
        std::random_device random_device;
        m_directory = std::filesystem::temp_directory_path() /
                      (prefix + "-" + std::to_string(unique) + "-" +
                       std::to_string(random_device()));
        std::filesystem::create_directories(m_directory);
    }

    ~TemporaryDirectory() {
        std::error_code error_code; // never throw from destructor
        std::filesystem::remove_all(m_directory, error_code);
    }

    TemporaryDirectory(const TemporaryDirectory &) = delete;
    TemporaryDirectory &operator=(const TemporaryDirectory &) = delete;

    /**
     * @brief writes content to a file in the temporary directory
     * @return path to the written file
     */
    std::filesystem::path write_file(const std::string &filename,
                                     const std::string &content) const {
        const auto filepath = m_directory / filename;
        std::ofstream file(filepath, std::ios::binary);
        assert(file.is_open());
        file << content;
        return filepath;
    }

    const std::filesystem::path &path() const { return m_directory; }

  private:
    std::filesystem::path m_directory;
};
