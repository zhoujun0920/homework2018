//
//  sudoku.h
//  homework2
//
//  Created by jzhou on 9/13/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#ifndef sudoku
#define sudoku
#include <vector>
#include <fstream>

using std::vector;
using namespace std;

class Sudoku {
    // Private
    int puzzle[9][9];
    int solution_num = 0;
    // Private member function that checks if the named row is valid
    bool row_valid(int row) {
        // write code that checks if "row" is valid
        for (int i = 0; i < 9; i++) {
            int temp = puzzle[row][i];
            if (temp != 0) {
                for (int j = i + 1; j < 9; j++) {
                    int k = puzzle[row][j];
                    if (k != 0 && temp == k) {
                        return false;
                    }
                }
            }
        }
        
        return true;
    }
    
    // Private member function that checks if the named column is valid
    
    bool col_valid(int col) {
        // check validity of "col"
        for (int i = 0; i < 9; i++) {
            int temp = puzzle[i][col];
            if (temp != 0) {
                for (int j = i + 1; j < 9; j++) {
                    int k = puzzle[j][col];
                    if (k != 0 && temp == puzzle[j][col]) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
    
    // Private member function that checks if the named 3x3 block is valid
    
    bool block_valid(int row, int col) {
        // check 3 x 3 block validity
        int tempRow = row / 3 * 3;
        int tempCol = col / 3 * 3;
        int temp = puzzle[row][col];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (tempRow + i != row && tempCol + j != col)
                    if (temp == puzzle[tempRow + i][tempCol + j]) {
                        return false;
                    }
            }
        }
        return true;
    }
    
public:
    // Public member function that reads the incomplete puzzle
    
    // we are not doing any checks on the input puzzle -- that is,
    
    // we are assuming they are indeed valid
    
    void read_puzzle(int argc, char* const argv[]) {
        // write code that reads the input puzzle using the
        
        // guidelines of figure 23 of the bootcamp material
        ifstream in(argv[1]);
        string fileName = argv[1];
        vector<char> p;
        char value_just_read_from_file;
        ifstream input_file(fileName);
        if (input_file) {
            while (input_file >> value_just_read_from_file) {
                p.push_back(value_just_read_from_file);
            }
            for (int i = 0; i < 9; i++) {
                for (int j = 0; j < 9; j++) {
                    puzzle[i][j] = p.at(i * 9 + j) - '0';
                }
            }
        } else {
            cout << "Couldn't open file" << endl;
        }
    }
    
    // Public member function that prints the puzzle when called
    
    void print_puzzle()
    
    {
        std::cout << "Board Position" << std::endl;
        
        for (int i = 0; i < 9; i++)
            
        {
            for (int j = 0; j < 9; j++)
                
            {
                // check if we have a legitimate integer between 1 and 9
                
                if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
                    
                {
                    // printing initial value of the puzzle with some formatting
                    
                    std::cout << puzzle[i][j] << " ";
                    
                }
                
                else {
                    // printing initial value of the puzzle with some formatting
                    
                    std::cout << "X ";
                }
            }
            
            std::cout << std::endl;
        }
    }
    
    // Public member function that (recursively) implements the brute-force
    
    // search for possible solutions to the incomplete Sudoku puzzle
    
    bool Solve(int row, int col) {
        // this part of the code identifies the row and col number of the
        
        // first incomplete (i.e. 0) entry in the puzzle.  If the puzzle has
        
        // no zeros, the variable row will be 9 => the puzzle is done, as
        
        // each entry is row-, col- and block-valid...
        
        // use the pseudo code of figure 3 of the description
        if (row == 9) return true;
        int nextRow = row, nextCol = col;
        if (col == 8) {
            nextRow++;
            nextCol = 0;
        } else {
            nextCol++;
        }
        if (puzzle[row][col] != 0) {
            return Solve(nextRow, nextCol);
        }
        for (int i = 1; i < 10; i++) {
            puzzle[row][col] = i;
            if (row_valid(row) && col_valid(col) && block_valid(row, col) &&
                Solve(nextRow, nextCol)) {
                solution_num++;
                cout << std::endl << "Solution #" << solution_num << endl;
                print_puzzle();
                continue;
            }
        }
        puzzle[row][col] = 0;
        return false;
    }
    
    void alternate_Solve(int row, int col) {
        
    }
};

#endif
