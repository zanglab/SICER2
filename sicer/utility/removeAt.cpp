/*
 * Code from: https://github.com/artem-ogre/remove_at/blob/master/remove_at.hpp
 * 
 * MIT License
 * 
 * Copyright (c) 2018 Artem Amirkhanov
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <algorithm>
#include <iterator>

/*!
 * Remove elements in the range [first; last) with indices from the sorted
 * unique range [ii_first, ii_last)
 */
template <class ForwardIt, class SortUniqIndsFwdIt>
inline ForwardIt remove_at(
    ForwardIt first,
    ForwardIt last,
    SortUniqIndsFwdIt ii_first,
    SortUniqIndsFwdIt ii_last)
{
    if(ii_first == ii_last) // no indices-to-remove are given
        return last;
    typedef typename std::iterator_traits<ForwardIt>::difference_type diff_t;
    typedef typename std::iterator_traits<SortUniqIndsFwdIt>::value_type ind_t;
    ForwardIt destination = first + static_cast<diff_t>(*ii_first);
    while(ii_first != ii_last)
    {
        // advance to an index after a chunk of elements-to-keep
        for(ind_t cur = *ii_first++; ii_first != ii_last; ++ii_first)
        {
            const ind_t nxt = *ii_first;
            if(nxt - cur > 1)
                break;
            cur = nxt;
        }
        // move the chunk of elements-to-keep to new destination
        const ForwardIt source_first =
            first + static_cast<diff_t>(*(ii_first - 1)) + 1;
        const ForwardIt source_last =
            ii_first != ii_last ? first + static_cast<diff_t>(*ii_first) : last;
        std::move(source_first, source_last, destination);
        // std::copy(source_first, source_last, destination) // c++98 version
        destination += source_last - source_first;
    }
    return destination;
}
