#pragma once

class FeatureToggle
{
    private:
        bool garbage_collect_;
    public:
        FeatureToggle(/* args */);
        ~FeatureToggle();

        FeatureToggle(bool garbage_collect)
        {
            garbage_collect_ = garbage_collect;
        }
        inline bool do_garbage_collection() const;
};

inline bool FeatureToggle::do_garbage_collection() const {
    return garbage_collect_;
}