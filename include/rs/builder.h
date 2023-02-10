#pragma once

#include <cassert>
#include <cmath>
#include <limits>

#include "common.h"
#include "radix_spline.h"

namespace rs {

// 构建器类
// 通过对有序数据的一次遍历来构建一个 RadixSpline 
// Allows building a `RadixSpline` in a single pass over sorted data.
template <class KeyType>
class Builder {
 public:
  Builder(KeyType min_key, KeyType max_key, size_t num_radix_bits = 18,
          size_t max_error = 32)
      : min_key_(min_key),
        max_key_(max_key),
        num_radix_bits_(num_radix_bits),
        num_shift_bits_(GetNumShiftBits(max_key - min_key, num_radix_bits)),
        max_error_(max_error),
        curr_num_keys_(0),
        curr_num_distinct_keys_(0),
        prev_key_(min_key),
        prev_position_(0),
        prev_prefix_(0) {
    // Initialize radix table, needs to contain all prefixes up to the largest
    // key + 1.
    // 见论文Page3: 要取r-bit的前缀,会有2^r个条目,要为radix_table预留 (2^r + 1)的空间, 该值=max_prefix+2;
    const uint32_t max_prefix = (max_key - min_key) >> num_shift_bits_;
    std::cout << "max_prefix=" << max_prefix << "; r=" << num_radix_bits_ 
              << "; shift=" << num_shift_bits_ << std::endl;
    radix_table_.resize(max_prefix + 2, 0);
  }

  // 遍历原始数据时逐个添加,但只有部分点会加入spline中
  // Adds a key. Assumes that keys are stored in a dense array.
  void AddKey(KeyType key) {
    // 如果当前key个数为0,就从开始添加;否则在前一个位置的后面添加
    if (curr_num_keys_ == 0) {
      AddKey(key, /*position=*/0);
      return;
    }
    AddKey(key, prev_position_ + 1);
  }

  // 完成构造并返回一个只读的 RS
  // Finalizes the construction and returns a read-only `RadixSpline`.
  RadixSpline<KeyType> Finalize() {
    // Last key needs to be equal to `max_key_`.
    assert(curr_num_keys_ == 0 || prev_key_ == max_key_);

    // Ensure that `prev_key_` (== `max_key_`) is last key on spline.
    if (curr_num_keys_ > 0 && spline_points_.back().x != prev_key_)
      AddKeyToSpline(prev_key_, prev_position_);

    // Maybe even size the radix based on max key right from the start
    FinalizeRadixTable();

    return RadixSpline<KeyType>(
        min_key_, max_key_, curr_num_keys_, num_radix_bits_, num_shift_bits_,
        max_error_, std::move(radix_table_), std::move(spline_points_));
  }

 private:
  // 获取需要移位的位数，基于最大key和最小key之间的差值 diff
  // KeyType == uint32_t 时
  // Returns the number of shift bits based on the `diff` between the largest
  // and the smallest key. KeyType == uint32_t.
  static size_t GetNumShiftBits(uint32_t diff, size_t num_radix_bits) {
    const uint32_t clz = __builtin_clz(diff);      // diff 二进制形式左起0的个数

    // 32-clz 就把 diff 二进制形式中前面的0全去掉了，剩下的就是其有效位数
    // 如 4 = 0b 0000 0000, 0000 0000, 0000 0000, 0000 0100; clz=29; 32-clz=3, 即4的有效位就是0b100中的后3位
    // 再让有效位数跟前缀 num_radix_bits比较, 若有效位数比要求的前缀还少,即只有3位,却要求前缀是5位,显然是不可能的,故设为0
    // 若有效位数超过了要求的前缀后,它俩的差值就是需要移位的位数；如：论文图1中的47183=0b 1011 1000 0100 1111, 32-clz=32-16=16; 且设置了前缀r=3，因此需要右移16-3=13位;
    if ((32 - clz) < num_radix_bits)
    {
      return 0;
    }
    return 32 - num_radix_bits - clz;
  }
  // KeyType == uint64_t.
  static size_t GetNumShiftBits(uint64_t diff, size_t num_radix_bits) {
    const uint32_t clzl = __builtin_clzl(diff);
    if ((64 - clzl) < num_radix_bits) return 0;
    return 64 - num_radix_bits - clzl;
  }

  // 在指定位置添加key
  void AddKey(KeyType key, size_t position) {
    assert(key >= min_key_ && key <= max_key_);
    // Keys need to be monotonically increasing. keys 需要单调递增
    assert(key >= prev_key_);
    // Positions need to be strictly monotonically increasing. positions 需要严格递增
    assert(position == 0 || position > prev_position_);

    PossiblyAddKeyToSpline(key, position);

    ++curr_num_keys_;
    prev_key_ = key;
    prev_position_ = position;
  }

  // 把点{key, position}添加到样条曲线上;就是放进 spline_points_数组中，以及把信息同步到 radix_table中
  // 把点添加到表示曲线的spline_points_数组中时，因为先更新该数组，再更新 radix_table，故在radix_table中key就直接映射到 spline_points_中最后一个元素；
  void AddKeyToSpline(KeyType key, double position) {
    spline_points_.push_back({key, position});
    PossiblyAddKeyToRadixTable(key);
  }

  // 方向：共线，顺时针，逆时针
  enum Orientation { Collinear, CW, CCW };
  static constexpr double precision = std::numeric_limits<double>::epsilon();   // 浮点数可表示的最小值（最小的浮点数）

  static Orientation ComputeOrientation(const double dx1, const double dy1,
                                        const double dx2, const double dy2) {
    const double expr = std::fma(dy1, dx2, -std::fma(dy2, dx1, 0));
    if (expr > precision)
      return Orientation::CW;
    else if (expr < -precision)
      return Orientation::CCW;
    return Orientation::Collinear;
  };

  void SetUpperLimit(KeyType key, double position) {
    upper_limit_ = {key, position};
  }
  void SetLowerLimit(KeyType key, double position) {
    lower_limit_ = {key, position};
  }

  // 设置前一个CDF点
  void RememberPreviousCDFPoint(KeyType key, double position) {
    prev_point_ = {key, position};
  }

  // Implementation is based on `GreedySplineCorridor` from:
  // T. Neumann and S. Michel. Smooth interpolating histograms with error
  // guarantees. [BNCOD'08]
  void PossiblyAddKeyToSpline(KeyType key, double position) {
    if (curr_num_keys_ == 0) {
      // Add first CDF point to spline. 向样条曲线中添加第一个CDF点
      AddKeyToSpline(key, position);
      ++curr_num_distinct_keys_;
      RememberPreviousCDFPoint(key, position);  // 添加结束时将当前{点,位置}设置为前一个CDF点
      return;
    }

    if (key == prev_key_) {
      // No new CDF point if the key didn't change.
      return;
    }

    // New CDF point.
    ++curr_num_distinct_keys_;

    // 用第2个CDF点初始化 upper_limit_ 和 lower_limit_
    if (curr_num_distinct_keys_ == 2) {
      // Initialize `upper_limit_` and `lower_limit_` using the second CDF
      // point.
      SetUpperLimit(key, position + max_error_);
      SetLowerLimit(key, (position < max_error_) ? 0 : position - max_error_);
      RememberPreviousCDFPoint(key, position);
      return;
    }

    // `B` in algorithm. 取spline_points最后一个点作为基准点Base point
    const Coord<KeyType>& last = spline_points_.back();

    // Compute current `upper_y` and `lower_y`.
    // 计算当前key（即点C-Current）的y值上下界限
    const double upper_y = position + max_error_;
    const double lower_y = (position < max_error_) ? 0 : position - max_error_;

    // Compute differences.
    assert(upper_limit_.x >= last.x);
    assert(lower_limit_.x >= last.x);
    assert(key >= last.x); 
    const double upper_limit_x_diff = upper_limit_.x - last.x;
    const double lower_limit_x_diff = lower_limit_.x - last.x;
    const double x_diff = key - last.x;

    assert(upper_limit_.y >= last.y);
    assert(position >= last.y);
    const double upper_limit_y_diff = upper_limit_.y - last.y;
    const double lower_limit_y_diff = lower_limit_.y - last.y;
    const double y_diff = position - last.y;

    // last是文献22图1中的B点,当前点(key,position)是C点,prev_point_是B和C中间那个点
    // `prev_point_` is the previous point on the CDF and the next candidate to
    // be added to the spline. Hence, it should be different from the `last`
    // point on the spline.
    assert(prev_point_.x != last.x);

    // 是否要收窄error corridor: 若当前点位于error corridor范围内,则收窄;否则将prev_point_添加到spline中
    // Do we cut the error corridor?
    if ((ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff,     /* 如果线段BC在BU左侧 */
                            y_diff) != Orientation::CW) ||
        (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff,     /* 或如果线段BC在BL右侧 */
                            y_diff) != Orientation::CCW)) {
      // 即线段BC不在 error corridor 范围内时,将CDF的前一个点(B和C中间的点)加入spline中,它即将作为新的基准点B
      // Add previous CDF point to spline. 
      AddKeyToSpline(prev_point_.x, prev_point_.y);

      // Update limits. 更新误差上下界
      SetUpperLimit(key, upper_y);
      SetLowerLimit(key, lower_y);
    } else { // 若线段BC在error corridor范围内,则不向spline中添加点,并看情况更新误差上/下界
      assert(upper_y >= last.y);
      const double upper_y_diff = upper_y - last.y;
      // 若BU在BU'左侧,则更新误差上界(误差上下界朝着缩小的方向变化)
      if (ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff,
                             upper_y_diff) == Orientation::CW) {
        SetUpperLimit(key, upper_y);
      }
      
      // 若BL在BL'右侧,则更新误差下界(因为新的误差下界更窄)
      const double lower_y_diff = lower_y - last.y;
      if (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff,
                             lower_y_diff) == Orientation::CCW) {
        SetLowerLimit(key, lower_y);
      }
    }

    RememberPreviousCDFPoint(key, position);
  }

  // 往 radix_table_ 表中添加key; possibly意为“可能”添加;因为如果当前key的前缀和前一个点的前缀一样时,就不会添加
  void PossiblyAddKeyToRadixTable(KeyType key) {
    // 当前key的前缀; key-min_key_应该表示相对于最小key的偏移后再移位
    const KeyType curr_prefix = (key - min_key_) >> num_shift_bits_;
    // 只有跟前一个点的前缀不一样才添加
    if (curr_prefix != prev_prefix_) {
      const uint32_t curr_index = spline_points_.size() - 1;    // 索引从下标0开始,因此最近一次的索引值就是 size()-1；curr_index表示当前key对应的是曲线上的哪个点

      // radix table里存储的是一个个的前缀值（有点类似于用bit位做标志）
      // 从 prev_prefix_+1到当前的前缀，统统赋值成该key在样条曲线中的位置索引（多个前缀都指向同一个点,多对一）
      for (KeyType prefix = prev_prefix_ + 1; prefix <= curr_prefix; ++prefix)
        radix_table_[prefix] = curr_index;
      prev_prefix_ = curr_prefix;
    }
  }

  void FinalizeRadixTable() {
    ++prev_prefix_;
    const uint32_t num_spline_points = spline_points_.size();
    for (; prev_prefix_ < radix_table_.size(); ++prev_prefix_)
      radix_table_[prev_prefix_] = num_spline_points;
  }

  const KeyType min_key_;
  const KeyType max_key_;
  const size_t num_radix_bits_;           // 前缀位数
  const size_t num_shift_bits_;           // 需要移动的位数
  const size_t max_error_;

  std::vector<uint32_t> radix_table_;
  std::vector<Coord<KeyType>> spline_points_;     // 样条曲线上的点（数组里是一个个的点 {key, position}）

  size_t curr_num_keys_;
  size_t curr_num_distinct_keys_;         // 当前不同keys的个数
  KeyType prev_key_;
  size_t prev_position_;
  KeyType prev_prefix_;

  // Current upper and lower limits on the error corridor of the spline.
  Coord<KeyType> upper_limit_;    // 误差上界; 误差上下界中间的空隙就组成了 error corridor
  Coord<KeyType> lower_limit_;    // 误差下界

  // Previous CDF point.
  Coord<KeyType> prev_point_;     // 记录前一个CDF点,原始数据中相同key的点不会记录进prev_point_中;如果原始数据前3个都是1,则该值只会记录第一个点(1,0)
};

}  // namespace rs