class NumMatrix:

    def __init__(self, rows, cols):
        self.mat = []
        for i in range(rows):
            self.mat.append([])
            for j in range(cols):
                self.mat[i].append(0.0)

    def __getitem__(self, n):
        return self.mat[n]
        
    def num_rows (self):
        return len(self.mat)
    
    def num_cols (self):
        return len(self.mat[0])
    
    def get_value (self, i, j):
        if i>j: return self.mat[i][j]
        else: return self.mat[j][i]
    
    def set_value(self, i, j, value):
        if i>j: self.mat[i][j] = value
        else: self.mat[j][i] = value
    
    def print_mat(self):
        for r in self.mat: print(r)
        print()
    
    def min_dist_indexes (self):
        m = self.mat[1][0]
        res= (1,0)
        for i in range(1,self.num_rows()):
            for j in range(i):
                if self.mat[i][j] < m:
                    m = self.mat[i][j]
                    res = (i, j)
        return res
    
    def add_row(self, newrow):
        self.mat.append(newrow)

    def add_col(self, newcol):
        for r in range(self.num_rows()):
            self.mat[r].append(newcol[r])

    def remove_row(self, ind):
        del self.mat[ind]

    def remove_col(self, ind):
        for r in range(self.num_rows()):
            del self.mat[r][ind]

    def copy(self):
        newm = NumMatrix(self.num_rows(), self.num_cols())
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                newm.mat[i][j] = self.mat[i][j]
        return newm
    