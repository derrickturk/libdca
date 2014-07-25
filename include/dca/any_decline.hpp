#ifndef ANY_DECLINE_HPP
#define ANY_DECLINE_HPP

#include <memory>
#include <type_traits>
#include <utility>
#ifdef __GNUC__
#ifdef __GXX_RTTI
#include <typeinfo>
#endif
#endif
#ifndef DCA_NO_IOSTREAMS
#include <iostream>
#endif

namespace dca {

class any {
    public:
        template<class Decline>
        any(Decline&& d);

        any(const any& other);
        any(any&& other) noexcept;

        any& operator=(const any& other);
        any& operator=(any&& other) noexcept;
        any& operator=(any& other) // block template operator=
        {
            return *this = const_cast<const any&>(other);
        }

        template<class Decline>
        any& operator=(Decline&& d);

#ifdef __GNUC__
#ifdef __GXX_RTTI
        const std::type_info& type() const;
#endif
#endif

        double rate(double time) const;
        double cumulative(double time) const;

#ifndef DCA_NO_IOSTREAMS
        friend std::ostream& operator<<(std::ostream& os, const any& d);
#endif

    private:
        struct any_impl_base {
            virtual any_impl_base* clone() const = 0;

#ifdef __GNUC__
#ifdef __GXX_RTTI
            virtual const std::type_info& type() const = 0;
#endif
#endif

            virtual double rate(double time) const = 0;
            virtual double cumulative(double time) const = 0;

#ifndef DCA_NO_IOSTREAMS
            virtual std::ostream& stream_to(std::ostream& os) const = 0;
#endif

            virtual ~any_impl_base() noexcept {}
        };

        template<class Decline>
        struct any_impl final : public any_impl_base {
            any_impl(const Decline& d) : d_(d) {}
            any_impl(Decline&& d) : d_(std::move(d)) {}

            any_impl* clone() const override
            {
                return new any_impl(d_);
            }

#ifdef __GNUC__
#ifdef __GXX_RTTI
            const std::type_info& type() const override
            {
                return typeid(Decline);
            }
#endif
#endif

            double rate(double time) const override
            {
                return d_.rate(time);
            }

            double cumulative(double time) const override
            {
                return d_.cumulative(time);
            }

#ifndef DCA_NO_IOSTREAMS
            std::ostream& stream_to(std::ostream& os) const override
            {
                return os << d_;
            }
#endif

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

#ifdef __GNUC__
#ifdef __GXX_RTTI
inline const std::type_info& any::type() const
{
    return impl_->type();
}
#endif
#endif

#ifndef DCA_NO_IOSTREAMS
inline std::ostream& operator<<(std::ostream& os, const any& d)
{
    return d.impl_->stream_to(os);
}
#endif

}

#endif
