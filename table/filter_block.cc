// Copyright (c) 2012 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "table/filter_block.h"

#include "leveldb/filter_policy.h"
#include "util/coding.h"

namespace leveldb {

// See doc/table_format.md for an explanation of the filter block format.

// Generate new filter every 2KB of data
// 每2KB数据构建一个filter的entry，对外还是一个data block一个filter
static const size_t kFilterBaseLg = 11;
static const size_t kFilterBase = 1 << kFilterBaseLg;

FilterBlockBuilder::FilterBlockBuilder(const FilterPolicy* policy)
    : policy_(policy) {}

void FilterBlockBuilder::StartBlock(uint64_t block_offset) {
  // 第一次调用 block_offset 为0，啥都不做
  uint64_t filter_index = (block_offset / kFilterBase);
  assert(filter_index >= filter_offsets_.size());
  while (filter_index > filter_offsets_.size()) {
    GenerateFilter();
  }
}

void FilterBlockBuilder::AddKey(const Slice& key) {
  Slice k = key;
  start_.push_back(keys_.size());
  keys_.append(k.data(), k.size());
}

Slice FilterBlockBuilder::Finish() {
  // 1. 写入剩余的key
  if (!start_.empty()) {
    GenerateFilter();
  }

  // Append array of per-filter offsets
  // 2. 最后要补充一些信息：
  //  filter内部每个的偏移量
  const uint32_t array_offset = result_.size();
  for (size_t i = 0; i < filter_offsets_.size(); i++) {
    PutFixed32(&result_, filter_offsets_[i]);
  }

  // filter总大小
  PutFixed32(&result_, array_offset);
  // 内部拆分的阈值
  result_.push_back(kFilterBaseLg);  // Save encoding parameter in result
  return Slice(result_);
}

void FilterBlockBuilder::GenerateFilter() {
  const size_t num_keys = start_.size();
  // 1. 当前没有等待记录的key，直接添加就行
  if (num_keys == 0) {
    // Fast path if there are no keys for this filter
    filter_offsets_.push_back(result_.size());
    return;
  }

  // Make list of keys from flattened key structure
  // 2. 恢复真实的key，放到tmp_keys_
  //   keys_存的是所有键值拼接起来的字符串，并没有分开存储。
  //   而是用start_来确定每个key的起始位置
  start_.push_back(keys_.size());  // Simplify length computation
  tmp_keys_.resize(num_keys);
  for (size_t i = 0; i < num_keys; i++) {
    // 这里获得第i个键在keys_中的起始位置指针
    const char* base = keys_.data() + start_[i];
    size_t length = start_[i + 1] - start_[i];
    tmp_keys_[i] = Slice(base, length);
  }

  // Generate filter for current set of keys and append to result_.
  // 3. 构建一个新的filter
  filter_offsets_.push_back(result_.size());
  policy_->CreateFilter(&tmp_keys_[0], static_cast<int>(num_keys), &result_);

  tmp_keys_.clear();
  keys_.clear();
  start_.clear();
}

FilterBlockReader::FilterBlockReader(const FilterPolicy* policy,
                                     const Slice& contents)
    : policy_(policy), data_(nullptr), offset_(nullptr), num_(0), base_lg_(0) {
  size_t n = contents.size();
  if (n < 5) return;  // 1 byte for base_lg_ and 4 for start of offset array
  // 按照格式，最后是分割阈值
  base_lg_ = contents[n - 1];
  // 读4字节，解析array_offset，也就是写入偏移量总大小
  uint32_t last_word = DecodeFixed32(contents.data() + n - 5);
  if (last_word > n - 5) return;
  // 每2k的所对应filter的偏移量和大小
  data_ = contents.data();
  offset_ = data_ + last_word;
  // 每个偏移数组元素占四个字节
  num_ = (n - 5 - last_word) / 4;
}

// 根据解析出来的数据，然后进行判断filter是否匹配
bool FilterBlockReader::KeyMayMatch(uint64_t block_offset, const Slice& key) {
  // block_offset / 2^base_lg_
  // 计算出当前在哪个filter中
  uint64_t index = block_offset >> base_lg_;
  if (index < num_) {
    // 读4字节，即filter offset
    uint32_t start = DecodeFixed32(offset_ + index * 4);
    uint32_t limit = DecodeFixed32(offset_ + index * 4 + 4);
    if (start <= limit && limit <= static_cast<size_t>(offset_ - data_)) {
      Slice filter = Slice(data_ + start, limit - start);
      return policy_->KeyMayMatch(key, filter);
    } else if (start == limit) {
      // Empty filters do not match any keys
      return false;
    }
  }
  return true;  // Errors are treated as potential matches
}

}  // namespace leveldb
