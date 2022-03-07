#ifndef __THREAD_H_
#define __THREAD_H_

#include <atomic>              // std::atomic
#include <chrono>              // std::chrono
#include <condition_variable>  // std::condition_variable
#include <functional>          // std::bind
#include <mutex>               // std::mutex
#include <queue>               // std::queue
#include <thread>              // std::thread, std::this_thread
#include <utility>             // std::forward
#include <vector>              // std::vector

class Task {
 public:
  template <typename Func_T, typename... ARGS>
  Task(Func_T func, ARGS... args) {
    this->func = std::bind(func, std::forward<ARGS>(args)...);
  }
  void run() {
    this->func();
    return;
  }
  std::function<void()> func;
};

class ThreadPool {
 public:
  ThreadPool(size_t n) {
    is_running = true;
    running_task_num = 0;
    for (size_t i = 0; i < n; i++) {
      threads.push_back(new std::thread(&ThreadPool::thread_worker, this));
    }
  }

  void submit(Task *t) {
    std::unique_lock<std::mutex> lock(m_mutex);
    tasks.push(t);
    m_cond.notify_one();
    return;
  }

  void wait() {
    while (true) {
      if (tasks.empty() && (running_task_num == 0)) {
        break;
      }
      std::this_thread::sleep_for(std::chrono::microseconds(1000));
    }
  }

  ~ThreadPool() {
    do {
      is_running = false;
      std::unique_lock<std::mutex> lock(m_mutex);
      m_cond.notify_all();
    } while (0);
    for (size_t i = 0; i < threads.size(); i++) {
      threads[i]->join();
      delete threads[i];
    }
    return;
  }

 private:
  std::vector<std::thread *> threads;
  std::queue<Task *> tasks;
  std::atomic<int> running_task_num;
  std::mutex m_mutex;
  std::condition_variable m_cond;
  bool is_running;

  void thread_worker() {
    while (is_running) {
      Task *t = this->get();
      if (t == nullptr) break;
      running_task_num += 1;
      t->run();
      running_task_num -= 1;
      delete t;
    }
    return;
  }

  Task *get() {
    std::unique_lock<std::mutex> lock(m_mutex);
    while (is_running && tasks.empty()) m_cond.wait(lock);
    Task *t = nullptr;
    if (is_running) {
      t = tasks.front();
      tasks.pop();
    }
    return t;
  }
};

#endif
