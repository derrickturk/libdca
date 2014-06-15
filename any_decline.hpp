#ifndef ANY_DECLINE_HPP

#include <memory>
#include <type_traits>
#include <utility>

namespace dca {

class any {
    public:
        template<class Decline>
        any(Decline&& d);

        any(const any& other);
        any(any&& other) noexcept;

        any& operator=(const any& other);
        any& operator=(any&& other) noexcept;
        any& operator=(any& other) = delete; // block template operator=

        template<class Decline>
        any& operator=(Decline&& d);

        double rate(double time) const;
        double cumulative(double time) const;

    private:
        struct any_impl_base {
            virtual any_impl_base* clone() const = 0;

            virtual double rate(double time) const = 0;
            virtual double cumulative(double time) const = 0;

            virtual ~any_impl_base() noexcept {}
        }

        template<class Decline>
        struct any_impl : public any_impl_base {
            any_impl(const Decline& d) : d_(d) {}
            any_impl(Decline&& d) : d_(std::move(d)) {}

            any_impl* clone() const override final
            {
                return new any_impl(d_);
            }

            double rate(double time) const override final
            {
                return d_.rate(time);
            }

            double cumulative(double time) const override final
            {
                return d_.cumulative(time);
            }

            Decline d_;
        };

        std::unique_ptr<any_impl_base> impl_;
};

template<class Decline>
inline any::any(Decline&& d)
    : impl_(new any_impl<typename std::decay<Decline>::type>(
                std::forward<Decline>(d))) { }

inline any::any(const any& other)
    : impl_(other.impl_->clone()) { }

inline any::any(any&& other) noexcept
    : impl_(std::move(other.impl_)) { }

inline any& any::operator=(const any& other)
{
    impl_ = std::unique_ptr<any_impl_base>(other.impl_->clone());
    return *this;
}

inline any& any::operator=(any&& other) noexcept
{
    impl_ = std::move(other.impl_);
    return *this;
}

template<class Decline>
inline any& any::operator=(Decline&& d)
{
    impl_ = std::unique_ptr<any_impl_base>(
            new any_impl<typename std::decay<Decline>::type>(
                std::forward<Decline>(d)));
    return *this;
}

inline double any::rate(double time) const
{
    return impl_->rate(time);
}

inline double any::cumulative(double time) const
{
    return impl_->cumulative(time);
}

}

#define ANY_DECLINE_HPP
#endif
