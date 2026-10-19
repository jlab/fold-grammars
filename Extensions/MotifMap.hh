#pragma once
#include<iostream>
#include <numeric>
#include <optional>
#include<set>
#include <string>
#include<unordered_map>
#include<string_view>
#include<utility>
#include "ali_t.hh"
#include "sequence.hh"
#include "subsequence.hh"
#include <rtlib/rna.hh>

inline std::vector<long unsigned int> get_gaps(const char *s, unsigned int i, unsigned int j){
    std::vector<long unsigned int> gap_vec;
    for (long unsigned int pos = i; pos < j; pos++){
        if (s[pos] == GAP_BASE){
            gap_vec.push_back(pos);
        }
    }
    return gap_vec;
}

struct std_set_hash {
    template<typename T1>
    size_t operator () (const std::set<T1> &input_set) const {
        std::size_t seed = input_set.size();
        for (const auto& iter: input_set) {
            seed ^= iter + 0x9e3779b9  + (seed << 6) + (seed >> 2); //NOLINT
        }
        return seed;
    }
};

struct Motif {std::string seq; char abb {'X'};};
enum class direction_type: std::uint8_t {forward,reverse,both};

class CounterMap{
    //Simple Unordered Map to count number of occurences of letters in a vector

    public:
    //using Counter = std::unordered_map<const char, unsigned int>;
    CounterMap() = default;

    CounterMap(const std::vector<char>  input_vector){
        for (unsigned int i = 0; i < input_vector.size(); i++){
            if (Counter.find(input_vector[i]) == Counter.end()){
                Counter[input_vector[i]] = 1;
            }
            else{
                Counter[input_vector[i]]++;
            }
        }

    };

    bool empty(){
        return Counter.empty();
    }

    // Source - https://stackoverflow.com/a/55961383
    // Posted by Alecto, modified by community. See post 'Timeline' for change history
    // Retrieved 2026-07-14, License - CC BY-SA 4.0
    const auto findMaxValuePair()
    {
        if (Counter.empty()) {
            throw std::runtime_error("CounterMap is empty. Cannot find max value pair.");
        }
        const auto compare = [](const auto &p1, const auto &p2){
            return p1.second < p2.second;
        };
        return *std::max_element(Counter.begin(), Counter.end(), compare);
    }

    private:
    std::unordered_map<char, unsigned int> Counter{};
};


class MotifMap{
    
    public:
    using DupeHashMap = std::unordered_map<Basic_Sequence<char, unsigned int>, std::set<char>, Hash_ali_array>; //Type for the Dupe Map recording all possible motifs for each sequence
    using MotifHashMap = std::unordered_map<Basic_Sequence<char, unsigned int>, char, Hash_ali_array>; //Type for the central Motif Map that contains motifs with their abbreviations
    using TranslationHashMap = std::unordered_map<std::set<char>, char,std_set_hash>; //Type for the translation map that connects 
    direction_type direction;
    std::string_view input;
    

    MotifMap() = default;

    //Most used constructor to create a MotifHash Map from an input string view and a direction for the motifs
    MotifMap(const std::string_view input_string, direction_type dir): direction(dir),input(input_string){ //Constructor definition
        Dupes = create_DupeMap(input, direction);
        Translations = create_translations(Dupes);
        Motifs = create_motif_hashmap(Dupes,Translations);
    };

    
    MotifHashMap::iterator find(const Basic_Subsequence<char, unsigned int> &input_subsequence){
        Basic_Sequence Motif {&input_subsequence.front(),input_subsequence.size()};
        return Motifs.find(Motif);
    };

    template<typename alphabet, typename pos_type>
    MotifHashMap::iterator find(const Basic_Subsequence<alphabet,pos_type> &internal_subsequence1, const Basic_Subsequence<alphabet,pos_type> &internal_subsequence2){
        Basic_Sequence Motif1 {&internal_subsequence1.front(),internal_subsequence1.size()};
        Basic_Sequence Motif2 {&internal_subsequence2.front(),internal_subsequence2.size()};
        Motif1.concat(Motif2.seq,Motif2.size());
        return Motifs.find(Motif1);
    };

    std::vector<char> find(const Basic_Subsequence<M_Char, unsigned int> &input_subsequence){
        std::vector<char> found;
        for (unsigned row = 0; row < rows(input_subsequence); row++){
            Basic_Sequence Motif;
            std::vector<long unsigned int> vec = get_gaps(input_subsequence.seq->row(row),input_subsequence.i,input_subsequence.j);
            Motif = Basic_Sequence {input_subsequence.seq->row(row), input_subsequence.i,input_subsequence.j, vec};
            if (auto search = Motifs.find(Motif); search != Motifs.end()) {
                for (unsigned int i = 0; i < Dupes[Motif].size();i++){
                    char mots = *next(Dupes[Motif].begin(),i);
                    found.push_back(mots);
                }
            }
        }
        return found;
    }

    std::vector<char> find(const Basic_Subsequence<M_Char, unsigned int> &input_subsequence1, const Basic_Subsequence<M_Char, unsigned int> &input_subsequence2){
        std::vector<char> found;
        for (unsigned int row = 0; row < rows(input_subsequence1); row++){
            std::vector<long unsigned int> vec1 = get_gaps(input_subsequence1.seq->row(row),input_subsequence1.i,input_subsequence1.j);
            std::vector<long unsigned int> vec2 = get_gaps(input_subsequence2.seq->row(row),input_subsequence2.i,input_subsequence2.j);
            Basic_Sequence Motif {input_subsequence1.seq->row(row), input_subsequence1.i,input_subsequence1.j,vec1};
            Basic_Sequence Motif2 {input_subsequence2.seq->row(row), input_subsequence2.i,input_subsequence2.j,vec2};
            Motif.concat(Motif2.seq,Motif2.size());
            if (auto search = Motifs.find(Motif);search != Motifs.end()) {
                for (unsigned int i = 0; i < Dupes[Motif].size();i++){
                    char mots = *next(Dupes[Motif].begin(),i);
                    found.push_back(mots);
                }
            }
        }
        return found;
    }  

    MotifHashMap::iterator find(const Basic_Sequence<char,unsigned int> &input_sequence){
        return Motifs.find(input_sequence);
    };

    template<typename alphabet, typename pos_type>
    MotifHashMap::iterator find(const Basic_Sequence<alphabet, pos_type> &input_sequence, unsigned int start_pos, unsigned int end_pos){
        const Basic_Subsequence Subseq{input_sequence, start_pos, end_pos};
        Basic_Sequence Motif {&Subseq.front(),Subseq.size()};
        return Motifs.find(Motif);
    };

    template<typename alphabet, typename pos_type>
    MotifHashMap::iterator find(Basic_Sequence<alphabet,pos_type> &internal_subsequence1, Basic_Sequence<alphabet,pos_type> &internal_subsequence2){
            internal_subsequence1.concat(internal_subsequence2.seq,internal_subsequence2.size());
        return Motifs.find(internal_subsequence1);
    };

    MotifHashMap::iterator begin() {return Motifs.begin();};
    MotifHashMap::iterator end() {return Motifs.end();};

    MotifHashMap::const_iterator begin() const {return Motifs.begin();};
    MotifHashMap::const_iterator end() const {return Motifs.end();};

    void add_motifs(std::string_view new_input){
        switch (direction) {
            case direction_type::forward:
                add_dupe_entries(new_input, false,Dupes);
                break;
            case direction_type::reverse:
                add_dupe_entries(new_input, true,Dupes);
                break;
            case direction_type::both:
                add_dupe_entries(new_input, false,Dupes);
                add_dupe_entries(new_input, true,Dupes);
                break;
        }
        Translations = create_translations(Dupes);
        Motifs = create_motif_hashmap(Dupes,Translations);
    };

    void print_duplicates(){
        std::vector<std::string> strings;
        for (const auto &keyvalue: Translations){
            std::vector<std::string> substring_vector;
            if (keyvalue.first.size() > 1){
            for (char mot: keyvalue.first){
                std::string substring(1,mot);
                substring_vector.emplace_back(substring);
            }
            std::string sub = std::accumulate(std::next(substring_vector.begin()), substring_vector.end(),substring_vector[0],[](const std::string& Existing_string, const std::string& Added_string){ return Existing_string +"/" + Added_string;});
            sub += " -> ";
            sub += keyvalue.second;
            strings.emplace_back(sub);
        }}
        if (strings.size() > 0) {
            std::string dupe_string = std::accumulate(std::next(strings.begin()), strings.end(),strings[0],[](const std::string& Existing_string, const std::string& Added_string){ return Existing_string +", " + Added_string;});
            std::cerr << "Attention, the same sequence appears for different motifs. Ambiguity cases: " << dupe_string << "\n";
        }
    };

    DupeHashMap Dupes;
    TranslationHashMap Translations;
    MotifHashMap Motifs;

    //Initializer functions
    static DupeHashMap create_DupeMap(const std::string_view input, direction_type direction){
        DupeHashMap DupeMap;
        switch (direction) {
            case direction_type::forward:
                add_dupe_entries(input, false, DupeMap);
                break;
            case direction_type::reverse:
                add_dupe_entries(input, true, DupeMap);
                break;
            case direction_type::both:
                add_dupe_entries( input, false, DupeMap);
                add_dupe_entries( input, true, DupeMap);
                break;
        }
        return DupeMap;
    };
    
    static TranslationHashMap create_translations(const DupeHashMap& DupeMap){
        TranslationHashMap translation_map;
        std::set<std::set<char>> dupe_sets;
        std::set<char> used_characters;
        for (const auto& keyvalue: DupeMap) {
            dupe_sets.insert(keyvalue.second);
        }
        for (const std::set<char>& set: dupe_sets){
            bool inserted = false;
            for (const char option: set){
                if (set.size() == 1 ){
                    used_characters.insert(option);
                    translation_map[set] = option;
                    inserted = true;
                    break;
                }
                if (used_characters.find(std::tolower(option,std::locale())) == used_characters.end()){
                    used_characters.insert(std::tolower(option,std::locale()));
                    translation_map[set] = std::tolower(option,std::locale());
                    inserted = true;
                    break;
                }
                continue;
            if (!inserted){
                throw std::runtime_error("One of the sets can not be used because all characters in it are used for other sets");
            }
            }
            //This is where some sort of logic would be nice to switch around the keys to see if any other combination is possible
        }
        return translation_map;
    };

    static MotifHashMap create_motif_hashmap(const DupeHashMap& DupeMap, const TranslationHashMap& TranslationMap){
        MotifHashMap MotiMap;
        for (const auto& keyvalue: DupeMap){
            MotiMap[keyvalue.first] = TranslationMap.at(keyvalue.second);
        }
        return MotiMap;
    };

    //General purpose in class functions
    //Function to add entries into DupeHashMap 
    static void add_dupe_entries(const std::string_view input_string_v, bool dir_bool, DupeHashMap& DupeMap){
        size_t pos_b = 0;
        size_t pos_e = 0;
        while ((pos_e = input_string_v.find('\n',pos_b)) != std::string::npos) {
            split_parse_line(pos_b, pos_e, input_string_v, dir_bool, DupeMap);
            pos_b = pos_e + 1 ;
        }
        if (pos_b <= input_string_v.size()){ //This catches the edge case of the last line not having a new line to use as a separator
            split_parse_line(pos_b, input_string_v.size()-pos_b, input_string_v,  dir_bool, DupeMap);
        }
    };
    

    static void split_parse_line(const size_t begin, const size_t end, const std::string_view view, bool dir_bool, DupeHashMap& DupeMap){
        std::string line;
        line = view.substr(begin,end-begin);
        if (const auto motif = parse_motif(line, dir_bool)){
            Basic_Sequence<char, unsigned int> basic_seq_motif {motif.value().seq.data(), static_cast<unsigned int>(motif.value().seq.size())};
            char_to_rnali(basic_seq_motif);
            if (auto dupe_sarch = DupeMap.find(basic_seq_motif); dupe_sarch == DupeMap.end()){ //Check if the sequence is already in the map
                DupeMap[basic_seq_motif] = std::set<char> {motif.value().abb};
            }
            else {
                DupeMap[basic_seq_motif].insert(motif.value().abb);
            }
        }
    };

    static std::optional<Motif> parse_motif (const std::string &input, bool rev) {
        constexpr char delim = ',';
        std::stringstream sinput (input);
        size_t elements {0};
        Motif motif{};
        std::string item;
        while (std::getline(sinput,item,delim)) {
            if (elements == 0){
                if (rev){
                    std::reverse(item.begin(),item.end());
                }
                motif.seq = item;
            }
            else if (elements == 1){
                motif.abb = item.front();
            }
            else{
                std::cerr << "Error more than one motif in one line" << "\n";
                return std::nullopt;
            }
            ++elements;
        }
        if (elements > 2){
            std::cerr << "Found " << elements << ", Expected " << 2 << "\n";
            return std::nullopt;
        }
        if (motif.abb == 'X'){
            return std::nullopt;
        }
        return motif;
    };
};