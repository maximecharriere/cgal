// Copyright (c) 2019  GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Necip Fazil Yildiran

#ifndef CGAL_POINTMATCHER_TRANSPOSE_PROPERTIES_H
#define CGAL_POINTMATCHER_TRANSPOSE_PROPERTIES_H

#if defined(CGAL_LINKED_WITH_POINTMATCHER) || defined(DOXYGEN_RUNNING)

#include <CGAL/pointmatcher/compute_registration_transformation.h>

namespace CGAL
{

  namespace pointmatcher
  {

    template <typename Scalar>
    using Matcher = typename PointMatcher<Scalar>::Matcher;

    namespace internal
    {

      template <typename Scalar, typename NamedParameters1>
      std::shared_ptr<Matcher<Scalar>>
      construct_matcher(const NamedParameters1 &np1)
      {
        typedef PointMatcher<Scalar> PM;

        using parameters::choose_parameter;
        using parameters::get_parameter;
        // using CGAL::pointmatcher::ICP_config;

        std::shared_ptr<Matcher<Scalar>> matcher;

        const ICP_config null_config{"_null_pm_config_in_cgal"};
        auto is_null_config = [&](const ICP_config &c) { return !c.name.compare(null_config.name); };

        // Matcher
        auto matcher_config = choose_parameter(get_parameter(np1, internal_np::matcher), null_config);
        if (!is_null_config(matcher_config))
        {
          try
          {
            matcher = PM::get().MatcherRegistrar.create(matcher_config.name, matcher_config.params);
          }
          catch (typename PointMatcherSupport::InvalidElement &error)
          {
            dump_invalid_point_matcher_config_exception_msg(error);
          }
        }

        return matcher;
      }

      template <class Scalar,
                class PointRange1,
                class PointRange2,
                class PointMap1,
                class PointMap2,
                class VectorMap1,
                class VectorMap2,
                class PropertyMap1,
                class PropertyMap2>
      bool transpose_properties(const PointRange1 &range1, const PointRange2 &range2,
                                PointMap1 point_map1, PointMap2 point_map2,
                                VectorMap1 vector_map1, VectorMap2 vector_map2,
                                PropertyMap1 property_map1, PropertyMap2 &property_map2,
                                std::shared_ptr<Matcher<Scalar>> matcher)
      {

        using PM = PointMatcher<Scalar>;
        using PM_cloud = typename PM::DataPoints;
        using PM_matrix = typename PM::Matrix;
        using PM_transform = typename PM::Transformation;
        using PM_transform_params = typename PM::TransformationParameters;

        // typedef typename boost::property_traits<PropertyMap>::value_type PropertyType;

        // ref_points: 1, points: 2
        std::size_t nb_ref_points = range1.size();
        std::size_t nb_points = range2.size();

        PM_matrix ref_points_pos_matrix = PM_matrix(4, nb_ref_points);
        PM_matrix ref_points_normal_matrix = PM_matrix(3, nb_ref_points);
        PM_matrix points_pos_matrix = PM_matrix(4, nb_points);
        PM_matrix points_normal_matrix = PM_matrix(3, nb_points);

        // In CGAL, point_set_1 is the reference while point_set_2 is the data

        // convert cgal points to pointmatcher points
        copy_cgal_points_to_pm_matrix<Scalar>(range1,
                                              point_map1,
                                              vector_map1,
                                              ref_points_pos_matrix,     // out
                                              ref_points_normal_matrix); // out

        copy_cgal_points_to_pm_matrix<Scalar>(range2,
                                              point_map2,
                                              vector_map2,
                                              points_pos_matrix,     // out
                                              points_normal_matrix); // out

        PM_cloud ref_cloud = construct_PM_cloud<PM_cloud>(ref_points_pos_matrix, ref_points_normal_matrix);
        PM_cloud cloud = construct_PM_cloud<PM_cloud>(points_pos_matrix, points_normal_matrix);

        matcher->init(ref_cloud);
        PM::Matches matches = matcher->findClosests(cloud);

        int idx = 0;
        for (const auto &p : range2)
        {
          // put(PM, key, newValue)
          put(property_map2, p, get(property_map1, *(range1.begin()+matches.ids(idx))));
          ++idx;
        }

        return true;
      }

    } // end of namespace internal

    template <class PointRange1, class PointRange2,
              class NamedParameters1, class NamedParameters2>
    bool transpose_properties(const PointRange1 &point_set_1, const PointRange2 &point_set_2,
                              const NamedParameters1 &np1, const NamedParameters2 &np2)
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;

      namespace PSP = CGAL::Point_set_processing_3;

      // Test property map types
      typedef typename CGAL::GetPointMap<PointRange1, NamedParameters1>::type PointMap1;
      typedef typename CGAL::GetPointMap<PointRange2, NamedParameters2>::type PointMap2;
      CGAL_static_assertion_msg((boost::is_same<typename boost::property_traits<PointMap1>::value_type,
                                                typename boost::property_traits<PointMap2>::value_type>::value),
                                "The point type of input ranges must be the same");

      typedef typename PSP::GetNormalMap<PointRange1, NamedParameters1>::type NormalMap1;
      typedef typename PSP::GetNormalMap<PointRange2, NamedParameters2>::type NormalMap2;
      CGAL_static_assertion_msg((boost::is_same<typename boost::property_traits<NormalMap1>::value_type,
                                                typename boost::property_traits<NormalMap2>::value_type>::value),
                                "The vector type of input ranges must be the same");


      typedef typename PSP::GetK<PointRange1, NamedParameters1>::Kernel Kernel;
      typedef typename Kernel::FT Scalar;

      // Extract property
      PointMap1 point_map1 = choose_parameter(get_parameter(np1, internal_np::point_map), PointMap1());
      NormalMap1 normal_map1 = choose_parameter(get_parameter(np1, internal_np::normal_map), NormalMap1());
      PointMap2 point_map2 = choose_parameter(get_parameter(np2, internal_np::point_map), PointMap2());
      NormalMap2 normal_map2 = choose_parameter(get_parameter(np2, internal_np::normal_map), NormalMap2());
      auto property_map1 = get_parameter(np1, internal_np::property_map);
      auto property_map2 = get_parameter(np2, internal_np::property_map);

      return internal::transpose_properties<Scalar>(point_set_1, point_set_2,
                                                    point_map1, point_map2,
                                                    normal_map1, normal_map2,
                                                    property_map1, property_map2,
                                                    internal::construct_matcher<Scalar>(np1));
    }

  }
} // end of namespace CGAL::pointmatcher

#endif // CGAL_LINKED_WITH_POINTMATCHER

#endif // CGAL_POINTMATCHER_TRANSPOSE_PROPERTIES_H
