#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "common.h"

namespace rs {

// 使用样条插值逼近累积分布函数（CDF）
// Approximates a cumulative distribution function (CDF) using spline
// interpolation.
template <class KeyType>
class RadixSpline {
 public:
  RadixSpline() = default;

  RadixSpline(KeyType min_key, KeyType max_key, size_t num_keys,
              size_t num_radix_bits, size_t num_shift_bits, size_t max_error,
              std::vector<uint32_t> radix_table,
              std::vector<rs::Coord<KeyType>> spline_points)
      : min_key_(min_key),
        max_key_(max_key),
        num_keys_(num_keys),
        num_radix_bits_(num_radix_bits),
        num_shift_bits_(num_shift_bits),
        max_error_(max_error),
        radix_table_(std::move(radix_table)),
        spline_points_(std::move(spline_points)) {}

  // 获取key的预估位置
  // Returns the estimated position of `key`.
  double GetEstimatedPosition(const KeyType key) const {
    // Truncate to data boundaries.
    if (key <= min_key_) return 0;                // 小于最小key，预估位置为0
    if (key >= max_key_) return num_keys_ - 1;    // 大于等于最大key，预估位置为keys数目数-1

    // 寻找 spline segment ，使得 key 位于 (spline[index-1], spline[index])范围内
    // 找key所在样条线段（key位于该样条段的末端）
    // Find spline segment with `key` ∈ (spline[index - 1], spline[index]].
    const size_t index = GetSplineSegment(key);
    const Coord<KeyType> down = spline_points_[index - 1];      // key所在样条段的起点
    const Coord<KeyType> up = spline_points_[index];            // key所在样条段的终点

    // 根据两个点down和up,可以求出过这两点的直线的方程(先求出斜率,再根据点斜式求出方程)
    // 计算斜率 k = (y1 - y2) / (x1 - x2);
    // Compute slope.
    const double x_diff = up.x - down.x;
    const double y_diff = up.y - down.y;
    const double slope = y_diff / x_diff;

    // 像论文中图2一样,真实的情况是CDF曲线，而通过插值（通过已知的、离散的数据点，在范围内推求新数据点的过程或方法）
    // 通过上面已经能求出样条线段的方程，而key又位于这个范围内,因此将key作为x代入，即可得 y-y0 = k(x-x0), 即y=k(x-x0)+y0;
    // (x0,y0)就是点down, key_diff就是(x-x0), 因此最终返回值就是通过插值预估的key的位置(key相当于一个x坐标，这里返回的相当于key的y坐标）
    // Interpolate.
    const double key_diff = key - down.x;
    return std::fma(key_diff, slope, down.y);
  }

  // 获取预估位置附近的搜索界限 [begin, end)
  // Returns a search bound [begin, end) around the estimated position.
  SearchBound GetSearchBound(const KeyType key) const {
    const size_t estimate = GetEstimatedPosition(key);
    // key的预估位置减去误差; 如论文中图2所示,estimate只是key在直线中的位置,而实际位置在CDF曲线上,两条线在y轴方向的距离就是最大误差;
    // 因此estimate-max_error就正好对应key在CDF上的位置(论文中图2画的有问题，真正的虚线应该垂直于x轴,不过代码是正确的，见https://github.com/learnedsystems/RadixSpline/issues/1
    const size_t begin = (estimate < max_error_) ? 0 : (estimate - max_error_);
    // `end` is exclusive. 
    // TODO: 为何是 estimate + max_error_+2 ?
    const size_t end = (estimate + max_error_ + 2 > num_keys_)
                           ? num_keys_
                           : (estimate + max_error_ + 2);
    return SearchBound{begin, end};
  }

  // 获取 rs的大小（占多少bytes）
  // Returns the size in bytes.
  size_t GetSize() const {
    std::cout << "sizeof(*this): " << sizeof(*this) << "; radix_table_.size: " << radix_table_.size()
            << "; spline_points_.size: " << spline_points_.size() << std::endl;

    return sizeof(*this) + radix_table_.size() * sizeof(uint32_t) +
           spline_points_.size() * sizeof(Coord<KeyType>);
  }

 private:
  // 获取样条点 spline point 的索引，该样条点标志着包含key的样条段的结束，其中 key 位于 (spline[index-1], spline[index])范围内
  // Returns the index of the spline point that marks the end of the spline
  // segment that contains the `key`: `key` ∈ (spline[index - 1], spline[index]]
  size_t GetSplineSegment(const KeyType key) const {
    // 使用 radix table 来缩小搜索范围
    // TODO: 取前缀为何要让 key-min_key_? 移位的数值是多少？
    // Narrow search range using radix table.
    const KeyType prefix = (key - min_key_) >> num_shift_bits_;     // 取key的前缀:
    assert(prefix + 1 < radix_table_.size());
    // 对应论文中图1
    // TODO: redix_table_中的值如何来的？看谁通过构造函数传进来的
    const uint32_t begin = radix_table_[prefix];
    const uint32_t end = radix_table_[prefix + 1];

    // 根据范围大小，选择线性查找或二分查找
    // TODO: 阈值为何是32？Tlx里搜索时好像也有相关阈值。
    if (end - begin < 32) {
      // Do linear search over narrowed range. // 范围较小时，在已缩小的范围内直接进行线性查找
      uint32_t current = begin;
      while (spline_points_[current].x < key) 
      {
        // TODO: 如果spline_points_中所有点都不等于key呢？那返回的也是第一个大于key的点
        ++current;
      }
      return current;
    }

    // 范围超过32 使用二分查找
    // 寻找大于等于key的第一个点的位置
    // Do binary search over narrowed range.
    const auto lb = std::lower_bound(
        spline_points_.begin() + begin, spline_points_.begin() + end, key,
        [](const Coord<KeyType>& coord, const KeyType key) {
          return coord.x < key;
        });
    // 利用distance计算出 spline_points_.begin() 和 lb 之间元素的数量，即所要找的索引
    return std::distance(spline_points_.begin(), lb);
  }

  KeyType min_key_;
  KeyType max_key_;
  size_t num_keys_;
  size_t num_radix_bits_;   // 前缀
  size_t num_shift_bits_;   // 需要移位的位数
  size_t max_error_;

  std::vector<uint32_t> radix_table_;
  std::vector<rs::Coord<KeyType>> spline_points_;   // 样条上的点

  template <typename>
  friend class Serializer;
};

}  // namespace rs