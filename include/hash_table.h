#pragma once

#ifdef DEBUG
#include <stdio.h>
#endif

#include <assert.h>
#include <stdlib.h>

#include "sys_utils.h"

/* a trivial hasher (meant for arithmetic types) */
template <typename K> struct DefaultHasher
{
        static constexpr K empty_key = ~static_cast<K>(0);
        size_t hash(K key) const { return static_cast<size_t>(key); }
        bool is_empty(K key) const { return (key == empty_key); }
        bool is_equal(K key1, K key2) const { return (key1 == key2); }
};

/* Uses OA and does _not_ implement deletion
 * Hasher must implement hash, is_empty,
 * is_equal and empty.		       */
template <typename K, typename V, typename H = DefaultHasher<K>>
struct HashTable
{
      public:
        /* Methods */
        HashTable(size_t expected_nkeys = 8, H hasher = H());
        ~HashTable();
        size_t size() const;
        void clear();
        void reserve(size_t expected_nkeys);
        V *get(K key) const;
        V *get_or_set(K key, V alt_val);
        void set_at(K key, V val);
        float load_factor() const;

      protected:
        /* Members */
        size_t _size;
        size_t _buckets;
        K *keys;
        V *vals;
        H hasher;
        /* Methods */
        void grow(size_t buckets);
        bool load_factor_ok() const;
};

/* Free function not linked to the class but useful in similar contexts */
template <typename K, typename H>
static inline size_t hash_lookup(K *keys, size_t buckets, H hasher, K key)
{
        assert(((buckets - 1) & buckets) == 0);
        size_t mask = buckets - 1;

        size_t bucket = hasher.hash(key) & mask;

        for (size_t probe = 0; probe < buckets; probe++)
        {

	if (hasher.is_empty(keys[bucket]) ||
	    hasher.is_equal(keys[bucket], key))
	{
	        return bucket;
	}
	/* quadratic probing */
	bucket = (bucket + probe + 1) & mask;
        }

        /* we should never reach this point */
        assert(false && "Table is full !\n");
        return 0;
}

template <typename K, typename V, typename H>
HashTable<K, V, H>::HashTable(size_t expected_keys, H hasher)
    : _size(0), _buckets(1), hasher(hasher)
{
        while (_buckets < (3 * expected_keys / 2))
        {
	_buckets *= 2;
        }

        keys = static_cast<K *>(malloc(_buckets * sizeof(K)));
        vals = static_cast<V *>(malloc(_buckets * sizeof(V)));
        assert(keys != nullptr && vals != nullptr);

        clear();
}

template <typename K, typename V, typename H> HashTable<K, V, H>::~HashTable()
{
        _buckets = 0;
        _size = 0;

        free(keys);
        free(vals);
}

template <typename K, typename V, typename H>
inline size_t HashTable<K, V, H>::size() const
{
        return (_size);
}

template <typename K, typename V, typename H> void HashTable<K, V, H>::clear()
{
        for (size_t i = 0; i < _buckets; ++i)
        {
	keys[i] = hasher.empty_key;
	assert(hasher.is_empty(keys[i]));
        }
        _size = 0;
}

template <typename K, typename V, typename H>
void HashTable<K, V, H>::reserve(size_t expected_keys)
{
        size_t buckets = 1;
        while (buckets < (3 * expected_keys / 2))
        {
	buckets *= 2;
        }

        grow(buckets);
}

template <typename K, typename V, typename H>
inline V *HashTable<K, V, H>::get(K key) const
{
        size_t bucket = hash_lookup(keys, _buckets, hasher, key);
        return hasher.is_empty(keys[bucket]) ? nullptr : &vals[bucket];
}

template <typename K, typename V, typename H>
inline V *HashTable<K, V, H>::get_or_set(K key, V alt_val)
{
        size_t bucket = hash_lookup(keys, _buckets, hasher, key);
        if (hasher.is_empty(keys[bucket]))
        {
	keys[bucket] = key;
	vals[bucket] = alt_val;
	_size++;
	if UNLIKELY (!load_factor_ok())
	{
	        grow(2 * _buckets);
	        assert(load_factor_ok());
	}

	return nullptr;
        }
        else
        {
	return &vals[bucket];
        }
}

template <typename K, typename V, typename H>
inline void HashTable<K, V, H>::set_at(K key, V val)
{
        size_t bucket = hash_lookup(keys, _buckets, hasher, key);
        vals[bucket] = val;
        if (hasher.is_empty(keys[bucket]))
        {
	keys[bucket] = key;
	_size++;
	if UNLIKELY (!load_factor_ok())
	{
	        grow(2 * _buckets);
	        assert(load_factor_ok());
	}
        }
}

template <typename K, typename V, typename H>
float HashTable<K, V, H>::load_factor() const
{
        return static_cast<float>(_size) / _buckets;
}

template <typename K, typename V, typename H>
void HashTable<K, V, H>::grow(size_t new_buckets)
{
        if (new_buckets <= _buckets)
	return;

#ifdef DEBUG
        printf("HashTable Grow to %zu!\n", new_buckets);
#endif

        assert((new_buckets & (new_buckets - 1)) == 0);

        K *newk = (K *)malloc(new_buckets * sizeof(*newk));
        V *newv = (V *)malloc(new_buckets * sizeof(*newv));

        for (size_t i = 0; i < new_buckets; ++i)
        {
	newk[i] = hasher.empty_key;
        }

        for (size_t probe = 0; probe < _buckets; ++probe)
        {

	const K key = keys[probe];

	if (hasher.is_empty(key))
	        continue;

	size_t new_idx = hash_lookup(newk, new_buckets, hasher, key);
	assert(hasher.is_empty(newk[new_idx]));

	newk[new_idx] = key;
	newv[new_idx] = vals[probe];
        }

        free(keys);
        free(vals);
        keys = newk;
        vals = newv;
        _buckets = new_buckets;
}

template <typename K, typename V, typename H>
inline bool HashTable<K, V, H>::load_factor_ok() const
{
        /* 66% load factor limit */
        return (_buckets > _size + _size / 2);
}
