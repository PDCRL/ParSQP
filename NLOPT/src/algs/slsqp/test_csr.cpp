#include<iostream>
#include<vector>

/*
void merge_update(int start_one, int end_one, int start_two, int end_two, std::vector<int> &row_ptr, std::vector<int> &col_idx, std::vector<double> &data, double sm, int rdx)
{

    //printf("Inside merge_update\n");    
    int ctr = start_two;
    int i = start_one;

    while (i < end_one && ctr < end_two)
    {   
        if (col_idx[ctr] == col_idx[i])
        {
            data[ctr] += sm * data[i];
            i++;
            ctr++;
        }
        
        else if (col_idx[ctr] < col_idx[i])
        {
            ctr++;
        }

        else if (col_idx[ctr] > col_idx[i])
        {
            col_idx.insert(col_idx.begin() + (ctr - 1), col_idx[i]);
            data.insert(data.begin() + (ctr - 1), sm * data[i]);
            
            for (int j = rdx+1; j < row_ptr.size(); j++)
            {
                row_ptr[j] += 1;
            }

            i++;
        }

    }

    int c = 0;
    while ( i < end_one)
    {
	col_idx.insert(col_idx.begin() + end_two + c, col_idx[i]);
	data.insert(data.begin() + end_two + c, sm * data[i]);

    	for (int j = rdx+1; j < row_ptr.size(); j++)
    	{
	    row_ptr[j] += 1;
    	}

    	i++;
	c++;
    }

}
*/

void merge_update(int start_one, int end_one, int start_two, int end_two, std::vector<int> &row_ptr, std::vector<int> &col_idx, std::vector<double> &data, double sm, int rdx, int *lpivot, double *up)
{
    printf("Inside Merge Update\n");
    int range = end_two-start_two;
    int counter = 0;
    for (int k = start_two; k < end_two; k++) {
	if (col_idx[k] > *lpivot-1){
		start_two = k;
		break;
	}
	counter++;
    }

    if((counter == range) || (col_idx[start_two-1] != *lpivot-1)) {
            col_idx.insert(col_idx.begin() + start_two, *lpivot-1);
            data.insert(data.begin() + start_two, sm * *up);
            
	    end_two++;

            for (int j = rdx; j < row_ptr.size(); j++)
            {
                row_ptr[j] += 1;
            }
    }
    else {
            data[start_two-1] += sm * *up;
    }

    int i1 = start_one;
    int i2 = start_two;

    while (i1 < end_one)
    {
	if (i2 < end_two) {
            if (col_idx[i1] == col_idx[i2])
            {
                data[i2] += sm * data[i1];
                i1++;
                i2++;
            }
            else if (col_idx[i1] > col_idx[i2])
            {
                i2++;
            }
	    else
	    {
                col_idx.insert(col_idx.begin() + i2, col_idx[i1]);
                data.insert(data.begin() + i2, sm * data[i1]);
            	end_two++;
                for (int j = rdx; j < row_ptr.size(); j++)
                {
                    row_ptr[j] += 1;
                }
                i1++;
	    }
	}
        else
        {
            col_idx.insert(col_idx.begin() + i2, col_idx[i1]);
            data.insert(data.begin() + i2, sm * data[i1]);
   	    end_two++;
            for (int j = rdx; j < row_ptr.size(); j++)
            {
                row_ptr[j] += 1;
            }
            i1++;
	    i2++;
        }
    }

}


void print_vector(std::vector<int> &v)
{
    for (auto &i: v)
        std::cout<<i<<" ";
}

int main()
{
    int sparse_matrix[4][4] = {{0, 1, 0, 1},
				{0, 1, 1, 0},
				{1, 0, 1, 1},
				{0, 1, 0, 1}};

    for (int i = 0; i<4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			std::cout<<sparse_matrix[i][j]<<" ";
		}
	std::cout<<"\n";
	}

    std::cout<<'\n';

    std::vector<int> row_ptr;
    std::vector<int> col_idx;
    std::vector<double> data;

    int matrix[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    row_ptr.resize(4 + 1, 0);
    
    for (int i = 0; i < 4; i++) 
    {    
    for (int j = 0; j < 4; j++) 
    {    
        if (sparse_matrix[i][j] != 0 || i == j)
        {    
            row_ptr[i+1] += 1;
            col_idx.push_back(j);
            data.push_back(sparse_matrix[i][j]);
        }    
    }    
    if (i + 1 < 4)
        row_ptr[i+2] = row_ptr[i+1];

    }    

    std::cout<<"Before Merge Update\nrow_ptr: ";
    print_vector(row_ptr); std::cout<<"\n";
    std::cout<<"col_idx: ";
    print_vector(col_idx); std::cout<<"\n";

    for (const auto &i : data)
        std::cout<<i<<" ";

    std::cout<<"\n\n";

    double n = -1.0;
    int lp = 1;

    merge_update(1, 3, 3, 5, row_ptr, col_idx, data, 2, 2, &lp, &n);

    std::cout<<"After Merge Update\nrow_ptr: ";
    print_vector(row_ptr); std::cout<<"\n";
    std::cout<<"col_idx: ";
    print_vector(col_idx); std::cout<<"\n";

    for (const auto &i : data)
        std::cout<<i<<" ";

    std::cout<<"\n\n";


for (int i = 0; i < row_ptr.size() - 1; i++) 
    {    
        int c_start = row_ptr[i];
        int c_end = row_ptr[i+1];

        for (int j = c_start; j < c_end; j++) 
        {    
            double val = data[j];
            int col = col_idx[j];
            matrix[i][col] = val; 

        }    

    }    

    for (int i = 0; i<4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			std::cout<<matrix[i][j]<<" ";
		}
	std::cout<<"\n";
	}


    return 0;
}
