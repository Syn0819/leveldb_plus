// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "db/log_writer.h"

#include <cstdint>

#include "leveldb/env.h"
#include "util/coding.h"
#include "util/crc32c.h"

namespace leveldb {
namespace log {

// 获取类型的校验码
static void InitTypeCrc(uint32_t* type_crc) {
  for (int i = 0; i <= kMaxRecordType; i++) {
    char t = static_cast<char>(i);
    type_crc[i] = crc32c::Value(&t, 1);
  }
}

Writer::Writer(WritableFile* dest) : dest_(dest), block_offset_(0) {
  InitTypeCrc(type_crc_);
}

Writer::Writer(WritableFile* dest, uint64_t dest_length)
    : dest_(dest), block_offset_(dest_length % kBlockSize) {
  InitTypeCrc(type_crc_);
}

Writer::~Writer() = default;

// 向Log文件中添加新增的键值对（Slice的形式）
// 日志结构：
/*
 *  <--------------header--------------->
 * |           ｜          ｜            ｜               ｜  
 *  crc (4byte) len (2byte) type (1byte)       content
*/
Status Writer::AddRecord(const Slice& slice) {
  const char* ptr = slice.data(); // 获取指向slice的实际数据
  size_t left = slice.size();     // 写入长度

  // Fragment the record if necessary and emit it.  Note that if slice
  // is empty, we still want to iterate once to emit a single
  // zero-length record
  // 意思是 写入内容会进行分片，但是如果写入数据为空，也会进行一次迭代
  // 以写入一个长度为0的记录
  Status s;
  // 这里begin标记这条记录是否为第一次写入，即如果一个记录跨越多个块
  // 
  bool begin = true;
  do {
    const int leftover = kBlockSize - block_offset_; // 当前块剩余空间
    assert(leftover >= 0);
    // 如果当前块剩下的空间连每个记录的header都放不下，则需要一个新的块，并将当前块剩余空间全部置0
    if (leftover < kHeaderSize) {
      // Switch to a new block
      if (leftover > 0) {
        // Fill the trailer (literal below relies on kHeaderSize being 7)
        static_assert(kHeaderSize == 7, "");
        // 为什么是\x00\x00\x00\x00\x00\x00？ 
        // 最多只有6bytes剩余
        dest_->Append(Slice("\x00\x00\x00\x00\x00\x00", leftover)); 
      }
      block_offset_ = 0;
    }

    // Invariant: we never leave < kHeaderSize bytes in a block.
    assert(kBlockSize - block_offset_ - kHeaderSize >= 0);

    // 当前块去掉一个header的长度剩下的可用空间
    const size_t avail = kBlockSize - block_offset_ - kHeaderSize;
    // 当前块能够写入的数据 取决于 剩余内容和块剩余空间之中较小的值
    const size_t fragment_length = (left < avail) ? left : avail;

    RecordType type;
    // 判断当前写入内容是否刚好写入到当前块
    const bool end = (left == fragment_length);
    // 通过begin和end字段组合判断header类型
    if (begin && end) {
      type = kFullType; // 当前写入，完整写了一个块
    } else if (begin) {
      type = kFirstType;  // 当前写入，是当前块的第一个写入
    } else if (end) {
      type = kLastType;   // 当前写入，是当前块的最后一个写入
    } else {
      type = kMiddleType; // 当前写入，是当前块的中间一个写入
    }

    s = EmitPhysicalRecord(type, ptr, fragment_length);
    ptr += fragment_length;
    left -= fragment_length;
    begin = false;
  } while (s.ok() && left > 0);
  return s;
}

// 函数功能：将日志写入磁盘
Status Writer::EmitPhysicalRecord(RecordType t, const char* ptr,
                                  size_t length) {
  // 保证写入长度不大于65535
  assert(length <= 0xffff);  // Must fit in two bytes
  // 确保记录写入长度不超过块大小
  assert(block_offset_ + kHeaderSize + length <= kBlockSize);

  // Format the header
  char buf[kHeaderSize];
  // little-endian
  buf[4] = static_cast<char>(length & 0xff);
  buf[5] = static_cast<char>(length >> 8);
  // type
  buf[6] = static_cast<char>(t);

  // Compute the crc of the record type and the payload.
  uint32_t crc = crc32c::Extend(type_crc_[t], ptr, length);
  crc = crc32c::Mask(crc);  // Adjust for storage
  EncodeFixed32(buf, crc);

  // Write the header and the payload
  Status s = dest_->Append(Slice(buf, kHeaderSize));
  if (s.ok()) {
    s = dest_->Append(Slice(ptr, length));
    if (s.ok()) {
      // WAL立即写盘
      s = dest_->Flush();
    }
  }
  block_offset_ += kHeaderSize + length;
  return s;
}

}  // namespace log
}  // namespace leveldb
