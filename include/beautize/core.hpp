// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace beautize {
constexpr double bohr_to_angstrom = 0.529177210903;
constexpr double hartree_to_ev = 27.211386245988;
constexpr double pi = 3.141592653589793238462643383279502884;
using Vector = std::vector<double>;
struct Error : std::runtime_error { using std::runtime_error::runtime_error; };
struct TextLine { std::string text, ending; };
struct Atom {
    int number{};
    bool frozen{};
    std::string label;
    std::size_t line{};
    std::array<std::pair<std::size_t, std::size_t>, 3> xyz_spans{};
};
struct Bond { int i{}, j{}; double order{}; };
struct Document {
    std::vector<TextLine> lines;
    std::vector<Atom> atoms;
    Vector positions; // Angstrom; flat [x1,y1,z1,x2,y2,z2,...]
    std::vector<Bond> bonds;
    std::vector<std::pair<std::size_t, std::string>> modredundant;
    std::string route;
    int charge{}, multiplicity{1};
    std::string render(const Vector& positions) const;
    Vector bond_matrix() const;
};
Document parse_gjf(const std::string& text);
std::string read_text(const std::string& path);
void write_atomic(const std::string& path, const std::string& contents, bool overwrite);
std::string trim(std::string s);
std::vector<std::string> words(const std::string& s);
int parse_int(const std::string& s, const std::string& context);
double parse_real(std::string s, const std::string& context);
std::string element_symbol(int z);
std::vector<int> parse_selection(const std::string& s, int nat);

enum class Kind { Bond = 2, Angle = 3, Dihedral = 4 };
struct Constraint {
    Kind kind{Kind::Bond};
    std::vector<int> atoms; // zero-based internally; canonical orientation
    double target{}; // Angstrom for B; radians for A and D
    std::string source;
};
struct ConstraintValue { double value{}; Vector gradient; };
ConstraintValue evaluate_coordinate(Kind kind, const std::vector<int>& atoms, const Vector& x);
double coordinate_error(Kind kind, double value, double target);
std::string constraint_name(const Constraint& c);
struct ConstraintSet {
    Vector reference;
    std::vector<bool> frozen;
    std::vector<Constraint> internals;
    double tolerance{1e-8}; // Angstrom (B), radians (A/D)
    explicit ConstraintSet(const Document& doc, bool ignore_gjf = false);
    void add(Kind kind, std::vector<int> ids, std::optional<double> target, const std::string& source);
    void add_modredundant(const std::string& text, const Document& doc, const std::string& source);
    void add_cli(Kind kind, const std::string& expression);
    bool empty() const;
    int frozen_count() const;
    double max_error(const Vector& x) const;
    void enforce_frozen(Vector& x) const;
    void restore(Vector& x, int max_iterations = 120) const;
    Vector project(const Vector& x, const Vector& v, int* rank = nullptr) const;
};
struct Evaluation { double energy{}; Vector gradient; }; // eV and eV/Angstrom
class Gfnff {
    void* handle_{};
    std::vector<int> numbers_;
public:
    static void set_threads(int count);
    Gfnff(const Document& doc, const Vector& x, int printlevel = 0);
    ~Gfnff();
    Gfnff(const Gfnff&) = delete;
    Gfnff& operator=(const Gfnff&) = delete;
    Evaluation operator()(const Vector& x);
    std::vector<int> graph() const;
};
struct OptimizeOptions {
    int max_steps{500}, memory{12}, threads{1};
    double fmax{0.05}, max_step{0.10};
};
struct StepInfo { int step{}; double energy{}, fmax{}, constraint_error{}, step_size{}; };
struct OptimizeResult {
    Vector positions;
    std::vector<StepInfo> history;
    bool converged{};
    std::string reason;
    int evaluations{};
};
OptimizeResult optimize(Vector x, const ConstraintSet& constraints,
    const std::function<Evaluation(const Vector&)>& evaluate,
    const OptimizeOptions& options,
    const std::function<void(const StepInfo&, const Vector&)>& observer = {});
double dot(const Vector& a, const Vector& b);
double max_atom_norm(const Vector& x);
std::string xyz_frame(const Document& doc, const Vector& x, const std::string& comment);
std::string make_report(const Document& doc, const ConstraintSet& constraints,
    const OptimizeResult& result, const OptimizeOptions& options,
    const std::string& input, const std::string& output, double initial_projection_displacement);
} // namespace beautize
