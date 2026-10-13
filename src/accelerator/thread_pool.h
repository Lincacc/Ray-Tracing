#pragma once

#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

namespace LinRT {
    // 所有工作线程共用这一把 mutex。lock_guard / unique_lock 离开作用域时自动解锁，不用手写 unlock。
    // 需要交给 condition_variable::wait 的地方用 unique_lock：wait 会先解锁，线程睡觉，被叫醒后再锁上。
    // 只是改一下队列或计数的地方用 lock_guard，它不能中途解锁。
    class thread_pool {
    public:
        explicit thread_pool(unsigned thread_count = std::thread::hardware_concurrency()) {
            if(thread_count == 0) {
                thread_count = 1;
            }

            workers.reserve(thread_count);
            for(unsigned i = 0; i < thread_count; ++i) {
                // 线程启动后进入 run()，在 cv 上睡觉，直到队列里有任务。
                workers.emplace_back([this] { run(); });
            }
        }

        thread_pool(const thread_pool&) = delete;
        thread_pool& operator=(const thread_pool&) = delete;

        ~thread_pool() {
            {
                std::lock_guard<std::mutex> lock(mutex);
                stop = true;
            }

            // 叫醒所有睡在 cv 上的线程。它们看到 stop 且队列空就退出 run()。
            cv.notify_all();
            for(auto& worker : workers) {
                // 阻塞等待，直到目标子线程彻底执行完毕并退出。
                worker.join();
            }
        }

        template<typename F>
        void submit(F&& job) {
            {
                std::lock_guard<std::mutex> lock(mutex);
                jobs.emplace(std::forward<F>(job));
            }
            // 解锁后再叫醒一个线程，让它可以立刻拿到锁去取任务。
            cv.notify_one();
        }

        void wait() {
            std::unique_lock<std::mutex> lock(mutex);
            // 队列空还不够：最后一块可能已被取走，线程仍在 job() 里。active 是正在执行的任务数。
            // done_cv.wait 睡觉前会放开这把锁，否则工作线程无法在结束时把 active 减 1。
            done_cv.wait(lock, [this] { return jobs.empty() && active == 0; });
        }

    private:
        void run() {
            for(;;) {
                std::function<void()> job;
                {
                    // unique_lock 才能交给 cv.wait。wait 期间会先解锁，线程睡觉，被叫醒后再锁上。
                    // 睡觉时必须放开 mutex，submit 才能把下一块放进队列。
                    std::unique_lock<std::mutex> lock(mutex);
        
                    cv.wait(lock, [this] { return stop || !jobs.empty(); });
                    if(stop && jobs.empty()) {
                        return;
                    }
        
                    job = std::move(jobs.front());
                    jobs.pop();
        
                    // 任务已取走但还没执行。wait() 靠 active 知道还有线程在干活。
                    ++active;
                }
        
                // 花括号结束，unique_lock 析构并解锁。渲染在锁外执行，其他线程才能同时领下一块。
                job();
                {
                    // 这里只改计数，不睡觉，lock_guard 离开作用域就会解锁。
                    std::lock_guard<std::mutex> lock(mutex);
                    --active;
                    if(jobs.empty() && active == 0) {
                        done_cv.notify_all();
                    }
                }
            }
        }

        std::vector<std::thread> workers;
        std::queue<std::function<void()>> jobs;
        std::mutex mutex;
        std::condition_variable cv;
        std::condition_variable done_cv;
        unsigned active = 0;
        bool stop = false;
    };
}