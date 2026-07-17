/*
 * energy_force.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: i-bird
 */

#ifndef EXAMPLE_VECTOR_4_REORDER_ENERGY_FORCE_HPP_
#define EXAMPLE_VECTOR_4_REORDER_ENERGY_FORCE_HPP_

constexpr int velocity = 0;
constexpr int force = 1;


template<typename CellList> void calc_forces(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double sigma12, double sigma6)
{
	// update the Cell-list
	vd.updateCellList(NN);

	// Get an iterator over particles
	auto it2 = vd.getDomainIterator();

	// For each particle p ...
	while (it2.isNext())
	{
		// ... get the particle p
		auto p = it2.get();

		// Get the position xp of the particle
		Point<3,double> xp = vd.getPos(p);

		// Reset the forice counter
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		// Get an iterator over the neighborhood particles of p
		auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

		// For each neighborhood particle ...
		while (Np.isNext())
		{
			// ... q
			auto q = Np.get();

			// if (p == q) skip this particle
			if (q == p.getKey())	{++Np; continue;};

			// Get the position of p
			Point<3,double> xq = vd.getPos(q);

			// Get the distance between p and q
			Point<3,double> r = xp - xq;

			// take the norm of this vector
			double rn = norm2(r);

			// Calculate the force, using pow is slower
			Point<3,double> f = 24.0*(2.0 *sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;

			// we sum the force produced by q on p
			vd.template getProp<force>(p)[0] += f.get(0);
			vd.template getProp<force>(p)[1] += f.get(1);
			vd.template getProp<force>(p)[2] += f.get(2);

			// Next neighborhood
			++Np;
		}

		// Next particle
		++it2;
	}
}


template<typename CellList> double calc_energy(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double sigma12, double sigma6)
{
	double E = 0.0;

	// Update the cell-list
	vd.updateCellList(NN);

	// Get the iterator
	auto it2 = vd.getDomainIterator();

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,double> xp = vd.getPos(p);

		// Reset the force
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		// Get an iterator over the neighborhood of the particle p
		auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			// if p == q skip this particle
			if (q == p.getKey())	{++Np; continue;};

			// Get position of the particle q
			Point<3,double> xq = vd.getPos(q);

			// take the normalized direction
			double rn = norm2(xp - xq);

			// potential energy (using pow is slower)
			E += 4.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) );

			// Next neighborhood
			++Np;
		}

		// Kinetic energy of the particle given by its actual speed
		E +=   (vd.template getProp<velocity>(p)[0]*vd.template getProp<velocity>(p)[0] +
				vd.template getProp<velocity>(p)[1]*vd.template getProp<velocity>(p)[1] +
				vd.template getProp<velocity>(p)[2]*vd.template getProp<velocity>(p)[2]) / 2;

		// Next Particle
		++it2;
	}

	return E;
}


#endif /* EXAMPLE_VECTOR_4_REORDER_ENERGY_FORCE_HPP_ */
