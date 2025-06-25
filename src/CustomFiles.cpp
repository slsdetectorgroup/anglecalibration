#include "CustomFiles.hpp"

namespace aare {

CustomMythenFile::CustomMythenFile(const ssize_t rows, const ssize_t cols,
                                   const std::string &mode) {
    m_rows = rows;
    m_cols = cols;
    m_mode = mode;
}

CustomMythenFile::~CustomMythenFile() { m_file.close(); }

Frame CustomMythenFile::read_frame() {
    auto f = Frame(m_rows, m_cols, m_dtype);
    uint32_t *frame_buffer = reinterpret_cast<uint32_t *>(f.data());
    uint32_t strip_index, photon_count;
    while (m_file >> strip_index >> photon_count) {
        *frame_buffer = photon_count;
        ++frame_buffer;
    }
    return f;
}

void CustomMythenFile::read_into(std::byte *image_buf) {
    uint32_t strip_index, photon_count;
    while (m_file >> strip_index >> photon_count) {
        std::memcpy(image_buf, &photon_count, sizeof(photon_count));
        image_buf += sizeof(photon_count);
    }
}

size_t CustomMythenFile::bytes_per_frame() {
    return m_rows * m_cols * m_dtype.bytes(); // TODO do i want m_counts?
}

size_t CustomMythenFile::pixels_per_frame() { return m_rows * m_cols; }

const Dtype CustomMythenFile::m_dtype(
    Dtype::TypeIndex::UINT32); // TODO: might not need or need it for python
                               // bindings

} // namespace aare