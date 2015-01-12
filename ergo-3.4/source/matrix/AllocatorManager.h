/* Ergo, version 3.4, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2014 Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */
#ifndef MAT_ALLOCATORMANAGER_HEADER
#define MAT_ALLOCATORMANAGER_HEADER

#include <stdexcept>
#include <list>
#include <string>
#include <sstream>
#include <iomanip> /* For setprecision */
#include "Allocator.h"

namespace mat {

template<class Treal>
class AllocatorManager
{
 public:
  void init(size_t noOfRealsPerBuffer_,
	    size_t noOfBuffers_) {
    if(noOfRealsPerBuffer != 0) {
      // This means that the AllocatorManager has already been initialized.
      // We allow this if the parameters are the same.
      if(noOfRealsPerBuffer_ != noOfRealsPerBuffer || noOfBuffers_ != noOfBuffers)
	throw std::runtime_error("Error in AllocatorManager: "
				 "attempt to re-initialize with different parameters.");
    }
    if(noOfRealsPerBuffer_ <= 0 || noOfBuffers_ <= 0)
      throw std::runtime_error("Error in AllocatorManager: bad input to init().");
    noOfRealsPerBuffer = noOfRealsPerBuffer_;
    noOfBuffers = noOfBuffers_;
  }
  static AllocatorManager & instance();
  Treal* alloc(size_t n) {
    if(n != noOfRealsPerBuffer)
      return new Treal[n];
    pthread_mutex_lock(&mutex);
    // Go through list to see if there is any free space.
    typename std::list< Allocator<Treal>* >::iterator it = list.begin();
    while(it != list.end()) {
      if(!(*it)->isFull()) {
	// OK, found allocator that is not full. Use it.
	Treal* ptr = (*it)->alloc();
	pthread_mutex_unlock(&mutex);
	return ptr;
      }
      it++;
    }
    // We did not find any non-full existing allocator. Need to add a new one.
    Allocator<Treal>* newAllocator = new Allocator<Treal>(noOfRealsPerBuffer, 
							 noOfBuffers);
    list.push_back(newAllocator);
    if(list.size() > peakListSize)
      peakListSize = list.size();
    Treal* ptr = newAllocator->alloc();
    pthread_mutex_unlock(&mutex);
    return ptr;
  }
  void free(Treal* ptr) {
    pthread_mutex_lock(&mutex);
    // Go through list to see if this ptr belongs to any allocator.
    typename std::list< Allocator<Treal>* >::iterator it = list.begin();
    while(it != list.end()) {
      if((*it)->ownsPtr(ptr)) {
	(*it)->free(ptr);
	// Now check if allocator is empty; in that case we want to remove it.
	if((*it)->isEmpty()) {
	  delete *it;
	  list.erase(it);
	}
	pthread_mutex_unlock(&mutex);
	return;
      }
      it++;
    }
    delete [] ptr;
    pthread_mutex_unlock(&mutex);
  }
  std::string getStatistics() {
    size_t noOfBytesPerAllocator = noOfBuffers * noOfRealsPerBuffer * sizeof(Treal);
    size_t totNoOfBytesAllocated = list.size() * noOfBytesPerAllocator;
    size_t peakNoOfBytesAllocated = peakListSize * noOfBytesPerAllocator;
    size_t totNoOfBytesUsed = 0;
    // Go through list to compute totNoOfBytesUsed
    typename std::list< Allocator<Treal>* >::iterator it = list.begin();
    while(it != list.end()) {
      totNoOfBytesUsed += (size_t)((*it)->getNoOfOccupiedSlots()) * noOfRealsPerBuffer * sizeof(Treal);
      it++;
    }
    std::stringstream ss;
    ss << "AllocatorManager statistics: ";
    ss << std::setprecision(3)
       << " noOfRealsPerBuffer: " << noOfRealsPerBuffer 
       << " noOfBuffers: " << noOfBuffers 
       << " list.size(): " << list.size() 
       << ". "
       << "Allocated: "  << (double)totNoOfBytesAllocated  / 1e9 << " GB, "
       << "Used: "       << (double)totNoOfBytesUsed       / 1e9 << " GB, "
       << "Peak alloc: " << (double)peakNoOfBytesAllocated/ 1e9 << " GB.";
    return ss.str();
  }
 private:
 AllocatorManager() : noOfRealsPerBuffer(0), noOfBuffers(0), peakListSize(0) {
    pthread_mutex_init(&mutex, NULL);
  }
  ~AllocatorManager() { 
    if(!list.empty())
      throw std::runtime_error("Error in AllocatorManager destructor: "
			       "not empty.");
    // Go through list to free any allocators that are left.
    typename std::list< Allocator<Treal>* >::iterator it = list.begin();
    while(it != list.end()) {
      delete *it;
      it++;
    }
  }
  std::list< Allocator<Treal>* > list;
  size_t noOfRealsPerBuffer;
  size_t noOfBuffers;
  pthread_mutex_t mutex;
  size_t peakListSize;
}; // end class AllocatorManager

} /* end namespace mat */

#endif
