#include "utils.hpp"
//--------------------------------------------------------------------------//
#include "globals.hpp"
//--------------------------------------------------------------------------//
using namespace std;
//--------------------------------------------------------------------------//
namespace dst::utils
{
  //--------------------------------------------------------------------------//
  void printSeparator()
  {
    cout << "---------------------------------------------------------------------" << endl;
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  using Clock_t = chrono::system_clock;
  using TimePoint_t = chrono::time_point<Clock_t>;
  using Duration_t = chrono::duration<double>;
  //--------------------------------------------------------------------------//
  string strTimePoint(const TimePoint_t &tp)
  {
    time_t time = Clock_t::to_time_t(tp);
    struct tm *p = localtime(&time);
    char raw[100];
    strftime(raw, sizeof(raw), "%d.%m.%Y, %H:%M:%S", p);
    return raw;
  }
  //--------------------------------------------------------------------------//
  void Timer::logCurrentTime(bool print_separator)
  {
    cout << "Current time: " << strTimePoint(Clock_t::now()) << endl;
    if (print_separator)
      printSeparator();
  }
  //--------------------------------------------------------------------------//
  Timer::Timer(bool start_active)
      : m_active(start_active),
        m_start_time(Clock_t::now()),
        m_prev_runtimes(Duration_t::zero()) {}
  //--------------------------------------------------------------------------//
  void Timer::start()
  {
    if (!m_active)
    {
      m_active = true;
      m_start_time = Clock_t::now();
    }
  }
  //--------------------------------------------------------------------------//
  void Timer::stop()
  {
    if (m_active)
    {
      m_active = false;
      m_prev_runtimes += Clock_t::now() - m_start_time;
    }
  }
  //--------------------------------------------------------------------------//
  void Timer::reset(bool set_active)
  {
    m_active = set_active;
    m_start_time = Clock_t::now();
    m_prev_runtimes = Duration_t::zero();
  }
  //--------------------------------------------------------------------------//
  const TimePoint_t &Timer::rawStartTime() const { return m_start_time; }
  //--------------------------------------------------------------------------//
  Duration_t Timer::rawRunTime() const
  {
    if (m_active)
      return Clock_t::now() - m_start_time + m_prev_runtimes;
    else
      return m_prev_runtimes;
  }
  //--------------------------------------------------------------------------//
  string Timer::strStartTime() const
  {
    return strTimePoint(m_start_time);
  }
  //--------------------------------------------------------------------------//
  string Timer::strRunTime() const
  {
    double total_seconds = rawRunTime().count();
    int minutes = (int)(total_seconds / 60.0);
    int seconds = (int)total_seconds % 60;
    int millis = (int)(total_seconds * 1000) % 1000;
    ostringstream stream;
    stream << minutes << "min " << seconds << "s " << millis << "ms";
    return stream.str();
  }
  //--------------------------------------------------------------------------//
  string Timer::strRunTimeCompact() const
  {
    double total_seconds = rawRunTime().count();
    int minutes = (int)(total_seconds / 60.0);
    int seconds = (int)total_seconds % 60;
    ostringstream stream;
    stream << '[' << (minutes < 10 ? "0" : "") << minutes
           << ':' << (seconds < 10 ? "0" : "") << seconds
           << ']';
    return stream.str();
  }
  //--------------------------------------------------------------------------//
  void Timer::logStartTime(bool print_separator) const
  {
    std::cout << "Start time: " << strStartTime() << endl;
    if (print_separator)
      printSeparator();
  }
  //--------------------------------------------------------------------------//
  void Timer::logRunTime(bool print_separator) const
  {
    std::cout << "Total time: " << strRunTime() << endl;
    if (print_separator)
      printSeparator();
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void RuntimeAnalyzer::startTimer(const string &name)
  {
    for (uint idx = 0; idx < m_timers.size(); ++idx)
    {
      if (m_timer_names[idx] == name)
      {
        m_timers[idx].start();
        return;
      }
    }
    throw runtime_error("RuntimeAnalyzer: no timer called '" + name + "'");
  }
  //--------------------------------------------------------------------------//
  void RuntimeAnalyzer::stopTimer(const string &name)
  {
    for (uint idx = 0; idx < m_timers.size(); ++idx)
    {
      if (m_timer_names[idx] == name)
      {
        m_timers[idx].stop();
        return;
      }
    }
    throw runtime_error("RuntimeAnalyzer: no timer called '" + name + "'");
  }
  //--------------------------------------------------------------------------//
  void RuntimeAnalyzer::logAnalysis(bool print_separator, RemnantHandling remnant_handling) const
  {
    cout << "Runtime Analyzer:" << '\n';
    if (m_timers.empty())
      cout << "(no timers)" << endl;
    else
    {
      // find maximum length of all names
      uint max_name_len = max_element(m_timer_names.begin(), m_timer_names.end(),
                                      [](const auto &s1, const auto &s2)
                                      { return s1.length() < s2.length(); })
                              ->length();
      // for PrintAsOther: min 5 required for 'Other'
      if (RemnantHandling::PrintAsOther == remnant_handling)
        max_name_len = max(5u, max_name_len);

      // sum up runtimes of all timers
      auto sumTimerTimes = [&timers = m_timers]() -> double
      {
        double sum_of_timers = 0.0;
        for (const Timer &t : timers)
          sum_of_timers += t.rawRunTime().count();
        return sum_of_timers;
      };

      // pretty print entries with their name and percentage
      auto printEntry = [max_name_len](const string &name, double time, double total_time) -> void
      {
        double percent = (total_time == 0.0 ? 0.0 : time / total_time * 100);
        stringstream stream;
        stream << fixed << setprecision(2) << setw(6) << percent << " %";
        prettyPrintListItem(name, stream.str(), max_name_len + 3);
      };

      // total time depending on RemnantHandling
      double total_time = RemnantHandling::Normalize == remnant_handling
                              ? sumTimerTimes()
                              : m_overall_timer.rawRunTime().count();
      // print all timings
      for (uint idx = 0; idx < m_timers.size(); ++idx)
      {
        double time = m_timers[idx].rawRunTime().count();
        printEntry(m_timer_names[idx], time, total_time);
      }
      if (RemnantHandling::PrintAsOther == remnant_handling)
        printEntry("Other", total_time - sumTimerTimes(), total_time);
    }

    if (print_separator)
      printSeparator();
  }
  //--------------------------------------------------------------------------//
  void RuntimeAnalyzer::reset()
  {
    for (Timer &timer : m_timers)
      timer.reset(false);
    m_overall_timer.reset(true);
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  LoopPrinter::LoopPrinter(const string &caption, uint iterations, bool end_with_separator)
      : m_caption(caption),
        m_current_iter(0),
        m_total_iters(iterations),
        m_end_with_separator(end_with_separator),
        m_finished_printing(false),
        m_timer(false)
  {
    // setup
    omp_init_lock(&m_lck);
    cout << caption << endl;
    m_timer.start();
    printProgress();
  }
  //--------------------------------------------------------------------------//
  LoopPrinter::~LoopPrinter()
  {
    omp_destroy_lock(&m_lck);
  }
  //--------------------------------------------------------------------------//
  void LoopPrinter::endIteration()
  {
    omp_set_lock(&m_lck);
    if (!m_finished_printing)
    {
      ++m_current_iter;
      printProgress();
    }
    omp_unset_lock(&m_lck);
  }
  //--------------------------------------------------------------------------//
  void LoopPrinter::printProgress()
  {
    cout << "\rSteps: " << m_current_iter << " / " << m_total_iters << flush;
    if (m_current_iter == m_total_iters)
    {
      cout << "\n"
           << m_timer.strRunTime() << endl;
      if (m_end_with_separator)
        printSeparator();
      m_finished_printing = true;
    }
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//