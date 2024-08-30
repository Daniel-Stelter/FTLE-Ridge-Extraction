#pragma once
//--------------------------------------------------------------------------//
#include "types.hpp"
//--------------------------------------------------------------------------//
#include <iomanip>
#include <omp.h>
#include <thread>
#include <vector>
//--------------------------------------------------------------------------//
namespace dst::utils
{
  //--------------------------------------------------------------------------//
  /** Helper function to suppress warnings due to unused parameters */
  template <class... Types>
  void unusedArgs(const Types... /*args*/) {}
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void printSeparator();
  //--------------------------------------------------------------------------//
  template <typename T>
  void prettyPrintListItem(const std::string &field, const T &val, uint field_space, bool line_break = true)
  {
    std::cout << "- " << std::setw(field_space) << std::left << field << val;
    if (line_break)
      std::cout << '\n';
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  class Timer
  {
    //--------------------------------------------------------------------------//
    using Clock_t = std::chrono::system_clock;
    using TimePoint_t = std::chrono::time_point<Clock_t>;
    using Duration_t = std::chrono::duration<double>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    static void logCurrentTime(bool print_separator = true);
    //--------------------------------------------------------------------------//
    Timer(bool start_active);
    //--------------------------------------------------------------------------//
    void start();
    void stop();
    void reset(bool set_active);
    //--------------------------------------------------------------------------//
    const TimePoint_t &rawStartTime() const;
    Duration_t rawRunTime() const;
    //--------------------------------------------------------------------------//
    std::string strStartTime() const;
    std::string strRunTime() const;
    std::string strRunTimeCompact() const;
    //--------------------------------------------------------------------------//
    void logStartTime(bool print_separator = true) const;
    void logRunTime(bool print_separator = true) const;
    //--------------------------------------------------------------------------//
  private:
    //--------------------------------------------------------------------------//
    bool m_active;
    TimePoint_t m_start_time;
    Duration_t m_prev_runtimes;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  class RuntimeAnalyzer
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    /**
     * Determines how the analysis should handle the margin between the total
     * execution time and the sum of the timers of all segments.
     */
    enum class RemnantHandling
    {
      Normalize,    // Percentages are normalized to add up to 100 %.
      PrintAsOther, // Remaining time is printed as additional segment 'Other'.
      Ignore        // Percentages may not add up to 100 %.
    };
    //--------------------------------------------------------------------------//
    template <typename... Args>
    RuntimeAnalyzer(Args... names)
        : m_overall_timer(true),
          m_timer_names{names...},
          m_timers(m_timer_names.size(), Timer(false)) {}
    //--------------------------------------------------------------------------//
    void startTimer(const std::string &name);
    //--------------------------------------------------------------------------//
    void stopTimer(const std::string &name);
    //--------------------------------------------------------------------------//
    void logAnalysis(bool print_separator = true, RemnantHandling remnant_handling = RemnantHandling::PrintAsOther) const;
    //--------------------------------------------------------------------------//
    void reset();
    //--------------------------------------------------------------------------//
  private:
    //--------------------------------------------------------------------------//
    Timer m_overall_timer;
    std::vector<std::string> m_timer_names;
    std::vector<Timer> m_timers;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  class LoopPrinter
  {
  public:
    //--------------------------------------------------------------------------//
    LoopPrinter(const std::string &caption,
                uint iterations,
                bool end_with_separator = true);
    //--------------------------------------------------------------------------//
    ~LoopPrinter();
    //--------------------------------------------------------------------------//
    void endIteration();
    //--------------------------------------------------------------------------//
  private:
    //--------------------------------------------------------------------------//
    void printProgress();
    //--------------------------------------------------------------------------//
    const std::string m_caption;
    uint m_current_iter;
    const uint m_total_iters;
    const bool m_end_with_separator;
    bool m_finished_printing;
    Timer m_timer;
    omp_lock_t m_lck;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//