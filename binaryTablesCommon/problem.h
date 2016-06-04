#ifndef PROBLEM_HEADER_GUARD
#define PROBLEM_HEADER_GUARD
#include <vector>
namespace binaryTables
{
	class problem
	{
	public:
		problem(std::vector<int>& rowSums, std::vector<int>& columnSums)
			:rowSums(rowSums), columnSums(columnSums)
		{}
		const std::vector<int>& getRowSums() const;
		const std::vector<int>& getColumnSums() const;
	private:
		std::vector<int> rowSums, columnSums;
	};
}
#endif
