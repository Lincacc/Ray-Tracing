#pragma once

#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

namespace LinRT {
    class thread_pool {
    public:
        explicit thread_pool(unsigned thread_count = std::thread::hardware_concurrency()) {
            if(thread_count == 0) {
                thread_count = 1;
            }

            workers.reserve(thread_count);
            for(unsigned i = 0; i < thread_count; ++i) {
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

            cv.notify_all();
            for(auto& worker : workers) {
                worker.join();
            }
        }

        template<typename F>
        void submit(F&& job) {
            {
                std::lock_guard<std::mutex> lock(mutex);
                jobs.emplace(std::forward<F>(job));
            }
            cv.notify_one();
        }

        void wait() {
            std::unique_lock<std::mutex> lock(mutex);
            done_cv.wait(lock, [this] { return jobs.empty() && active == 0; });
        }

    private:
        void run() {
            for(;;) {
                std::function<void()> job;
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    
                    cv.wait(lock, [this] { return stop || !jobs.empty(); });
                    if(stop && jobs.empty()) {
                        return;
                    }

                    job = std::move(jobs.front());
                    jobs.pop();

                    ++active;
                }

                job();
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    --active;
                    if(jobs.empty() && active == 0) {
                        done_cv.notify_all();
                    }
                }
            }
        }
    private:
        std::vector<std::thread> workers;
        std::queue<std::function<void()>> jobs;
        std::mutex mutex;
        std::condition_variable cv;
        std::condition_variable done_cv;
        unsigned active = 0;
        bool stop = false;
    };
}