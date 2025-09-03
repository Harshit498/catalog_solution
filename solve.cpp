#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

// Single-header JSON library by nlohmann
#include "json.hpp"

// Single-header BigInt library by matt-42
#include "BigInt.hpp"

// Use nlohmann's json and the BigInt library
using json = nlohmann::json;

// Function to convert a string representation of a number from any base to BigInt
BigInt stringToBigInt(const std::string& s, int base) {
    BigInt res = 0;
    BigInt p = 1;

    for (int i = s.length() - 1; i >= 0; i--) {
        int digit;
        if (s[i] >= '0' && s[i] <= '9') {
            digit = s[i] - '0';
        } else if (s[i] >= 'a' && s[i] <= 'f') {
            digit = s[i] - 'a' + 10;
        } else {
            throw std::invalid_argument("Invalid character in number string");
        }
        if (digit >= base) {
            throw std::invalid_argument("Digit exceeds base");
        }
        res = res + p * digit;
        p = p * base;
    }
    return res;
}

int main() {
    // Read the entire input into a JSON object
    json input;
    try {
        std::cin >> input;
    } catch (json::parse_error& e) {
        std::cerr << "JSON parsing error: " << e.what() << std::endl;
        return 1;
    }

    // Extract k, the number of points needed
    int k = input["keys"]["k"];

    // Store the points (x, y) as BigInts
    std::vector<std::pair<BigInt, BigInt>> points;
    
    // The keys "1", "2", etc. in the JSON are the x-values.
    // The inner 'base' and 'value' are used to calculate the y-values.
    for (auto& item : input.items()) {
        if (item.key() != "keys") {
            BigInt x = std::stoi(item.key());
            std::string val_str = item.value()["value"];
            int base = std::stoi(item.value()["base"].get<std::string>());
            BigInt y = stringToBigInt(val_str, base);
            points.push_back({x, y});
        }
    }
    
    // We only need k points to solve for the k coefficients
    if (points.size() < k) {
        std::cerr << "Error: Not enough points provided to solve the polynomial." << std::endl;
        return 1;
    }

    // Construct the augmented matrix [A|b] for the system Ax = b
    // where x is the vector of polynomial coefficients [a0, a1, ..., a(k-1)]
    std::vector<std::vector<BigInt>> matrix(k, std::vector<BigInt>(k + 1));

    for (int i = 0; i < k; ++i) {
        BigInt x = points[i].first;
        BigInt y = points[i].second;
        
        // The last column is the y value
        matrix[i][k] = y;
        
        // Fill the powers of x: x^0, x^1, x^2, ...
        BigInt x_pow = 1;
        for (int j = 0; j < k; ++j) {
            matrix[i][j] = x_pow;
            x_pow *= x;
        }
    }

    // --- Gaussian Elimination ---
    // Goal: Convert the matrix to row echelon form (upper triangular)
    for (int i = 0; i < k; ++i) {
        // Find pivot (non-zero element)
        int pivot_row = i;
        while (pivot_row < k && matrix[pivot_row][i] == 0) {
            pivot_row++;
        }

        if (pivot_row == k) {
            std::cerr << "Matrix is singular, cannot solve uniquely." << std::endl;
            return 1;
        }
        std::swap(matrix[i], matrix[pivot_row]);

        // Eliminate column i for all rows below row i
        for (int j = i + 1; j < k; ++j) {
            // To avoid fractions, we use cross-multiplication:
            // Row[j] = Row[j] * pivot_val - Row[i] * current_val
            BigInt pivot_val = matrix[i][i];
            BigInt current_val = matrix[j][i];
            
            for (int l = i; l <= k; ++l) {
                matrix[j][l] = matrix[j][l] * pivot_val - matrix[i][l] * current_val;
            }
        }
    }

    // --- Back Substitution ---
    // Goal: Solve for the coefficients from the bottom up
    std::vector<BigInt> coeffs(k);
    for (int i = k - 1; i >= 0; --i) {
        BigInt sum = 0;
        for (int j = i + 1; j < k; ++j) {
            sum += matrix[i][j] * coeffs[j];
        }
        
        // Solve for the coefficient a_i
        // matrix[i][i] * a_i = matrix[i][k] - sum
        coeffs[i] = (matrix[i][k] - sum) / matrix[i][i];
    }

    // The constant term is the first coefficient, a0
    std::cout << "The constant term (secret) is: " << coeffs[0] << std::endl;

    return 0;
}