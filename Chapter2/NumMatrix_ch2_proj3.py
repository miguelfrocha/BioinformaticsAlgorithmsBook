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
    

    def sum_mat (self):
        s = 0.0;
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                s +=  self.mat[i][j]        
        return s

    def mean (self):
        return self.sum_mat() / (self.num_rows() * self.num_cols())

    def maximum (self):
        m = self.mat[0][0]
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                if (self.mat[i][j] > m):
                    m = self.mat[i][j]        
        return m
    
    def minimum (self):
        m = self.mat[0][0]
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                if (self.mat[i][j] < m):
                    m = self.mat[i][j]        
        return m
    
    def sum_row (self, li):
        return sum(self.mat[li])
    
    def mean_row (self, li):
        s = 0
        for k in range(len(self.mat[li])):
            s += self.mat[li][k]
        return s / len(self.mat[li])

    def max_row (self, li):
        m = self.mat[li][0]
        for k in range(1,len(self.mat[li])):
            if self.mat[li][k] > m:
                m = self.mat[li][k]
        return m

    def min_row (self, li):
        m = self.mat[li][0]
        for k in range(1,len(self.mat[li])):
            if self.mat[li][k] < m:
                m = self.mat[li][k]
        return m

    def mean_rows(self):
        res = []
        for r in range(self.num_rows()):
            res.append(self.mean_row(r))
        return res

    def min_rows(self):
        res = []
        for r in range(self.num_rows()):
            res.append(self.min_row(r))
        return res

    def max_rows(self):
        res = []
        for r in range(self.num_rows()):
            res.append(self.max_row(r))
        return res

    def sum_col(self, lc):
        s = 0
        for k in range(self.num_rows()):
            s += self.mat[k][lc]
        return s

    def mean_col(self, lc):
        s = 0
        for k in range(self.num_rows()):
            s += self.mat[k][lc]
        return s / self.num_rows()
    
    def max_col(self, lc):
        m = self.mat[0][lc]
        for c in range(1, self.num_rows()):
            if self.mat[c][lc] > m:
                m = self.mat[c][lc]
        return m
    
    def min_col(self, lc):
        m = self.mat[0][lc]
        for c in range(1, self.num_rows()):
            if self.mat[c][lc] < m:
                m = self.mat[c][lc]
        return m
    
    def sum_cols(self):    
        res = []
        for c in range(self.num_cols()):
            res.append(self.sum_col(c))
        return res
    
    def mean_cols(self):
        res = []
        for c in range(self.num_cols()):
            res.append(self.mean_col(c))
        return res

    def max_cols(self):
        res = []
        for c in range(self.num_cols()):
            res.append(self.max_col(c))
        return res
    
    def min_cols(self):
        res = []
        for c in range(self.num_cols()):
            res.append(self.min_col(c))
        return res
    
    def square_mat(self):
        if self.num_cols() == self.num_rows():
            return True
        else: return False
    
    def mult_diagonal(self):
        if not self.square_mat(): 
            return None
        m = 1
        for r in range(self.num_rows()):
            m *= self.mat[r][r]
        return m
    
    def mult_scalar(self, scalar):
        m = NumMatrix(self.num_rows(), self.num_cols())
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                m.set_value(i, j, self.get_value(i,j) * scalar)
        return m
    
    def add_matrix(self, mat2):
        if mat2.num_rows() != self.num_rows() or mat2.num_cols() != self.num_cols():
            return None
        m = NumMatrix(self.num_rows(), self.num_cols())
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                m.set_value(i, j, self.get_value(i,j) + mat2.get_value(i, j) )
        return m       
    
    def mult_mat(self, mat2):
        if self.num_cols() != mat2.num_rows():
            return None
        m = NumMatrix(self.num_rows(), mat2.num_cols())
        for i in range(self.num_rows()):
            for j in range(mat2.num_cols()):
                s = 0
                for k in range(self.num_cols()):
                    s += self.get_value(i,k) * mat2.get_value(k,j)
                m.set_value(i,j,s)
        return m
    
def test():
    m = NumMatrix(3,3)
    m.set_value(0,0,1.0)
    m.set_value(2,1,3.0)
    m.set_value(1,0,2.0)
    m.set_value(0,1,5.0)
    m.set_value(1,1,-2.0)
    m.set_value(2,2,4.0)
    m.print_mat()
    print(m.mean())
    print(m.maximum())
    print("Mean rows")
    print(m.mean_rows())
    print("Mean cols")
    print(m.mean_cols())
    print("")
    print(m.square_mat())
    print(m.mult_diagonal())
    
    m1 = m.mult_scalar(5.0)
    m1.print_mat()
    
    m2 = m.add_matrix(m1)
    m2.print_mat()
    
    ident = NumMatrix(3,3)
    ident.set_value(0,0,1.0)
    ident.set_value(1,1,1.0)
    ident.set_value(2,2,1.0)
    ident.print_mat()
    
    m3 = m.mult_mat(ident)
    m3.print_mat()

if __name__ == '__main__':
    test()